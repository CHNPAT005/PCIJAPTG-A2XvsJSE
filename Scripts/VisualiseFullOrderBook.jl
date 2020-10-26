## Author: Patrick Chang & Ivan Jericevich
# Script file to visualise the full market depth on A2X

## Preamble

using CSV, CodecBzip2, DataFrames, ProgressMeter, Dates, Plots, Statistics, LaTeXStrings
cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")
# Create a dictionary mapping the securityIds to the security names
clearconsole()
SecurityIDtoTickerName = CSV.read("Supporting information/SecurityIDtoTickerName.csv")
secIDtoTickerName = Dict(SecurityIDtoTickerName[i, 1] => SecurityIDtoTickerName[i, 2] for i in 1:(size(SecurityIDtoTickerName)[1]))

## Phase 1 - Convert raw data to TAQ data
function getFromString(field::String, message::SubString{String}, delimiter::String) # Function to extract value of a specific field from a message string
    first_ind = findfirst(field*":", message)[end]+1
    last_ind = findfirst(delimiter, message[first_ind:end])[1]-2
    return message[first_ind:first_ind+last_ind]
end
function Unix2DateTime(time::Float64) # Function to convert unix time with nanoseconds into DateTime with milliseconds
    seconds = time ÷ 10^9
    milliseconds = Millisecond(time ÷ 10^6 - seconds * 10^3)
    nanoseconds = Nanosecond(time % 10^6)
    return unix2datetime(seconds + 2*3600) + milliseconds + nanoseconds
end
function makeDataStruct(message::SubString{String}, type::Int) # Function to create dataframe information from a message string. The structure of fields to enter depend on the type of message string.
    if type == 2
        securityId = getFromString("securityId", message, ",")
        side = getFromString("side", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("limitPrice", message, ","))
        orderId = getFromString("orderId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OA", price, quantity, side, "", "")
    elseif type == 3
        securityId = getFromString("securityId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        orderId = getFromString("orderId", message, ",")
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OC", NaN, NaN, "", "", "")
    elseif type == 4
        securityId = getFromString("securityId", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        orderId = getFromString("orderId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OM", price, quantity, "", "", "")
    elseif type == 5
        securityId = getFromString("securityId", message, ",")
        tradeType = getFromString("tradeType", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        orderId = getFromString("orderId", message, ",")
        tradeRef = getFromString("tradeRef", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "AT", price, quantity, "", tradeRef, tradeType)
    else # type 6 message
        securityId = getFromString("securityId", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        tradeRef = getFromString("tradeRef", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, "", "TB", price, NaN, "", tradeRef, "")
    end
end
function getTAQ(file::String) # Function to read in a day of data of all assets, remove useless information and construct dataframe of useful information for each asset. Each DataFrame is stored as a value in a dictionary.
    # Load in the data and split into vector
    rawData = read(file) |> x -> transcode(Bzip2Decompressor, x) |> y -> String(y) |> z -> split(z, "\n")
    # Filter only the relevant messages
    filter!(line -> occursin("msgType:2", line) || occursin("msgType:3", line) || occursin("msgType:4", line) || occursin("msgType:5", line) || occursin("msgType:6", line), rawData)
    # Truncate useless information
    messages = map(i -> split(rawData[i], " ")[2], 1:length(rawData))
    # Get the securityIds and names of the tickers active during this day
    securityIds = unique(getFromString.("securityId", messages, ","))
    tickerNames = map(i -> secIDtoTickerName["securityId:" * securityIds[i]], 1:length(securityIds))
    # Data Structure for each row of DataFrame of all assets
    taqAllSecurities = DataFrame([String, Float64, DateTime, String, String, Float64, Float64, String, String, String], [:Security, :TimeStamp, :Date, :OrderRef, :Type, :Price, :Quantity, :Side, :TradeRef, :TradeType])
    # Loop through the vector of messages
    for line in messages
        # Identify msgType
        msgType = parse(Int, getFromString("msgType", line, ","))
        # Extract TAQ format and push into dataframe
        getDataStruct = makeDataStruct(line, msgType)
        push!(taqAllSecurities, getDataStruct)
    end
    # Seperate securities by creating dictionary of DataFrames for each asset
    taqDict = Dict()
    for id in securityIds
        idIndeces = findall(x -> x == id, taqAllSecurities[:, :Security])
        if length(idIndeces) == 0
            emptyDataFrame = DataFrame([String, Float64, DateTime, String, String, Float64, Float64, String, String, String], [:Security, :TimeStamp, :Date, :OrderRef, :Type, :Price, :Quantity, :Side, :TradeRef, :TradeType])
            push!(taqDict, secIDtoTickerName[id] => emptyDataFrame)
        else
            push!(taqDict, secIDtoTickerName["securityId:" * id] => taqAllSecurities[idIndeces, :])
        end
    end
    return taqDict, securityIds, tickerNames
end

filedir = "Real Data/A2X/Raw/2019-07-15_A_prod_10G_stream.log.bz2"

NPN_CleanedData = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")
NPN_messageData = getTAQ(filedir)[1]["NPN"]


# Function to plot the Full Market depth for a given ticker
# Plots the MicroPrice, bid, ask, and trade over a user specified time of day.
# Requires both raw data and cleaned data.
# uses the clean data to plot trades and micro-price and raw data to build the
# order books for plotting.
function plotFullDepth(MessageData::DataFrame, CleanedData::DataFrame; date = "2019-07-15", start = 9, close = 17, startmin = 0, closemin = 0, pos = :topright)
    # Filter out data for the day
    dates = Date.(CleanedData[:,2])
    date = Date(date, "y-m-d")
    # Pull out the data for the day
    inds = findall(x -> x == date, dates)
    CleanedData = CleanedData[inds,:]


    # Get range of time to investigate
    StartTime = DateTime(Date(CleanedData[:Date][1])) + Hour(start) + Minute(startmin)
    EndTime = DateTime(Date(CleanedData[:Date][1])) + Hour(close) + Minute(closemin)
    # Get range of index to plot
    Inds = findall(x -> StartTime <= x <= EndTime, CleanedData[:Date])
    # Obtain range of data that matters
    CleanedData = CleanedData[Inds,:]


    # Get average Bid and Ask vol
    bidind = findall(!isnan, CleanedData[:BidVol])
    askind = findall(!isnan, CleanedData[:AskVol])
    tradeind = findall(!isnan, CleanedData[:TradeVol])

    # AveBidVol = mean(CleanedData[bidind, :BidVol])
    # AveAskVol = mean(CleanedData[askind, :AskVol])
    AveTradeVol = mean(CleanedData[tradeind, :TradeVol])
    AveBidVol = mean(MessageData[findall(x-> x == "BUY", MessageData[:Side]), :Quantity])
    AveAskVol = mean(MessageData[findall(x-> x == "SELL", MessageData[:Side]), :Quantity])



    # Plot MicroPrice
    p = plot(Time.(CleanedData[:Date]), CleanedData[:MicroPrice] .÷ 10^7, linetype = :steppost, color = :black, label = "Microprice", legend = pos, dpi = 300, linewidth = 0.5, size = (700, 400))

    ## Initialise dictionary for all bids and asks
    Bid_dict = Dict()
    Ask_dict = Dict()

    # Run through the TAQ messages to create dataframe of L1 BAT
    @showprogress "Computing..." for i in 1:size(MessageData)[1]
        line = MessageData[i,:]
        if line[:Date] < EndTime
            if line[:Type] == "OA"  # check if message is an Order Add
                if line[:Side] == "SELL"    # Order is ask
                    # Push data
                    push!(Ask_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Ask_dict)
                            # Update the plots
                            # Ask
                            plot!(p, [Time.(line[:Date])], [Ask_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :pink, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                else    # Order is bid
                    # Push data
                    push!(Bid_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Bid_dict)
                            # Update the plots
                            # Bid
                            plot!(p, [Time.(line[:Date])], [Bid_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :lightblue, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                end
            elseif line[:Type] == "OC"  # check if messeage is an Order Cancel
                # Check if the Order Cancel is in bid or ask
                if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask
                    # Update order book
                    delete!(Ask_dict, line[:OrderRef])
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Ask_dict)
                            # Update the plots
                            # Ask
                            plot!(p, [Time.(line[:Date])], [Ask_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :pink, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                else    # Order is from Bid
                    # Update order book
                    delete!(Bid_dict, line[:OrderRef])
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Bid_dict)
                            # Update the plots
                            # Bid
                            plot!(p, [Time.(line[:Date])], [Bid_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :lightblue, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                end
            elseif line[:Type] == "OM"  # check if messeage is an Order Modify
                # Check if the Order Modify is in bid or ask
                if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask
                    # Update order book
                    Ask_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Ask_dict)
                            # Update the plots
                            # Ask
                            plot!(p, [Time.(line[:Date])], [Ask_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :pink, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                else    # order is from Bid
                    # Update order book
                    Bid_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Only plot if it is in right time
                    if line[:Date] >= StartTime && line[:Date] <= EndTime
                        for j in keys(Bid_dict)
                            # Update the plots
                            # Bid
                            plot!(p, [Time.(line[:Date])], [Bid_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :lightblue, markerstrokealpha = 0.1, label = "", markerstrokewidth=0)
                        end
                    end
                end
            elseif line[:Type] == "AT"  # check if messeage is a Trade
                # Check if the Order Reference for trade is in bid or ask
                if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask; trade is buyer-initiated
                    # Minus trade volume from order volume
                    vol = Ask_dict[line[:OrderRef]][2] - line[:Quantity]
                    # Check if trade finished the order
                    if vol != 0 # trade did not finish entire order
                        # modify quantity in order
                        Ask_dict[line[:OrderRef]] = (Ask_dict[line[:OrderRef]][1], vol)
                        # Only plot if it is in right time
                        if line[:Date] >= StartTime && line[:Date] <= EndTime
                            for j in keys(Ask_dict)
                                # Update the plots
                                # Ask
                                plot!(p, [Time.(line[:Date])], [Ask_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :pink, markerstrokealpha = 0.01, label = "", markerstrokewidth=0)
                            end
                        end
                    else
                        # remove the order as it is finished
                        delete!(Ask_dict, line[:OrderRef])
                        # Only plot if it is in right time
                        if line[:Date] >= StartTime && line[:Date] <= EndTime
                            for j in keys(Ask_dict)
                                # Update the plots
                                # Ask
                                plot!(p, [Time.(line[:Date])], [Ask_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :pink, markerstrokealpha = 0.01, label = "", markerstrokewidth=0)
                            end
                        end
                    end
                elseif haskey(Bid_dict, line[:OrderRef]) # order is from Bid; trade is seller-initiated
                    # Minus trade volume from order volume
                    vol = Bid_dict[line[:OrderRef]][2] - line[:Quantity]
                    # Check if trade finished the order
                    if vol != 0 # trade did not finish entire order
                        # modify quantity in order
                        Bid_dict[line[:OrderRef]] = (Bid_dict[line[:OrderRef]][1], vol)
                        # Only plot if it is in right time
                        if line[:Date] >= StartTime && line[:Date] <= EndTime
                            for j in keys(Bid_dict)
                                # Update the plots
                                # Bid
                                plot!(p, [Time.(line[:Date])], [Bid_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :lightblue, markerstrokealpha = 0.01, label = "", markerstrokewidth=0)
                            end
                        end
                    else
                        # remove the order as it is finished
                        delete!(Bid_dict, line[:OrderRef])
                        # Only plot if it is in right time
                        if line[:Date] >= StartTime && line[:Date] <= EndTime
                            for j in keys(Bid_dict)
                                # Update the plots
                                # Bid
                                plot!(p, [Time.(line[:Date])], [Bid_dict[j][1]] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :lightblue, markerstrokealpha = 0.01, label = "", markerstrokewidth=0)
                            end
                        end
                    end
                end
            end
        end
    end
    current()

    # Make top-of-book blue and red
    # Bid
    plot!(p, Time.(CleanedData[bidind, :Date]), CleanedData[bidind, :Bid] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :blue, markerstrokealpha = 0.1, label = "Bid", markerstrokewidth=0)
    # Ask
    plot!(p, Time.(CleanedData[askind, :Date]), CleanedData[askind, :Ask] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :red, markerstrokealpha = 0.1, label = "Ask", markerstrokewidth=0)
    # Plot the trades separately
    plot!(p, Time.(CleanedData[tradeind, :Date]), CleanedData[tradeind, :Trade] .÷ 10^7, markersize = 1.5, seriestype = :scatter, color = :green, markerstrokealpha = 0.1, label = "Trade", markerstrokewidth=0)

    # Axis information
    f = "\\textrm{"
    xlabel!(p, latexstring("$f$(Date(CleanedData[:Date][1]))}"))
    ylabel!(p, L"\textrm{Price [ZAR]}")
    return p
end

mydate = "2019-07-15"

NPNFMDafternoon = plotFullDepth(NPN_messageData, NPN_CleanedData; date = mydate, start = 14, close = 14, startmin = 30, closemin = 40, pos = :outertopright)
# savefig(NPNFMDafternoon, "Plots/NPNFMDafternoon.png")

NPNFMDmorning = plotFullDepth(NPN_messageData, NPN_CleanedData; date = mydate, start = 13, close = 14, startmin = 50, closemin = 0, pos = :outertopright)
# savefig(NPNFMDmorning, "Plots/NPNFMDmorning.png")
