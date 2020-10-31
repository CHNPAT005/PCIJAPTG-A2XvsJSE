### Title: Data Cleaning - A2X
### Authors: Patrick Chang and Ivan Jericevich
### Function: Convert the raw A2X data into a detailed and usable L1LOB format.
### Structure:
# 1. Preliminaries
# 2. Cleaning functions
# 3. Implement cleaning functions
# 4. Calculate the frequencies of aggressive trades
### Strategy:
# Phase 1: Convert the raw message strings into a usable DataFrame format.
# Phase 2: Convert the cleaned messages into a usable L1 BAT order book with separate columns for each item.
# Phase 3: Create a more detailed dataframe with the addition of mid-price micro-price and inter-arrivals.
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, CodecBzip2, DataFrames, ProgressMeter, Dates
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
# Create a dictionary mapping the securityIds to the security names
SecurityIDtoTickerName = CSV.read("Supporting information/SecurityIDtoTickerName.csv")
secIDtoTickerName = Dict(SecurityIDtoTickerName[i, 1] => SecurityIDtoTickerName[i, 2] for i in 1:(size(SecurityIDtoTickerName)[1]))
#---------------------------------------------------------------------------


### 2. Cleaning Functions
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
function makeDataStruct(message::SubString{String}, type::Int) # Function to create dataframe information from a message string. The structure of fields to enter depend on the type of message string
    # Depending on the event type, the messages will have different fields
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
    else # Type 6 message
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
function combineTAQ(files) # Function to streamline everything. Reads in all the data and returns a dictionary of dataframes for the entire history of the A2X equities.
    MasterDictionary = Dict()
    @showprogress "Phase 1:" for file in files
        # Clean a single day of raw messages
        cleanDay = getTAQ(file); taq = cleanDay[1]; tickers = cleanDay[3]
        # Concactenate the resulting dictionary to the master dictionary
        for ticker in tickers
            if haskey(MasterDictionary, ticker)
                MasterDictionary[ticker] = vcat(MasterDictionary[ticker], taq[ticker])
            else
                push!(MasterDictionary, ticker => taq[ticker])
            end
        end
    end
    return MasterDictionary
end

## Phase 2 - Build L1 Order Book
function MakeL1BAT(data::DataFrame) # Function takes in 1 day of cleaned TAQ messages for 1 asset and creates the L1 Bid Ask Trade with their associated volumes. Trades are also classified into buyer- or seller-initiated.
    ## Initialise dictionary for all bids and asks
    Bid_dict = Dict(); Ask_dict = Dict()
    # Initialise DataFrame for L1 Order Book
    df = DataFrame(TimeStamp = Float64[], Date = DateTime[], EventType = String[], Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[], Trade = Float64[], TradeVol = Float64[], TradeSign = String[])
    # Run through the TAQ messages to create dataframe of L1 BAT
    for i in 1:size(data)[1]
        line = data[i,:]
        if line[:Type] == "OA"  # Check if message is an Order Add
            if line[:Side] == "SELL"
                # Check if Dictionary is empty
                if isempty(Ask_dict) # If empty simply push new ask into dictionary and df
                    push!(Ask_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, line[:Price], line[:Quantity], NaN, NaN, "")
                    push!(df, temp)
                else
                    # Get current best ask
                    best_ask = findmin(Ask_dict)[1][1]
                    # Push data
                    push!(Ask_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    # Get new ask
                    new_ask = line[:Price]
                    # Compare new ask against best ask
                    if new_ask < best_ask # If new ask is new best, update DataFrame
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, line[:Price], line[:Quantity], NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            else
                # Check if Dictionary is empty
                if isempty(Bid_dict) # If empty simply push new ask into dictionary and df
                    push!(Bid_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    temp = (line[:TimeStamp], line[:Date], "BID", line[:Price], line[:Quantity], NaN, NaN, NaN, NaN, "")
                    push!(df, temp)
                else
                    # Get current best bid
                    best_bid = findmax(Bid_dict)[1][1]
                    # Push data
                    push!(Bid_dict, line[:OrderRef] => (line[:Price], line[:Quantity]))
                    # Get new bid
                    new_bid = line[:Price]
                    # Compare new bid against best bid
                    if new_bid > best_bid # If new bid is new best, update DataFrame
                        temp = (line[:TimeStamp], line[:Date], "BID", line[:Price], line[:Quantity], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            end
        elseif line[:Type] == "OC" # Check if messeage is an Order Cancel
            # Check if the Order Cancel is in bid or ask
            if haskey(Ask_dict, line[:OrderRef]) # If true, order is from Ask
                # Get orderRef of best ask
                ref = findmin(Ask_dict)[2]
                if ref != line[:OrderRef] # If order to be removed is NOT best ask then remove order from dictionary. Don't need to print new best ask
                    delete!(Ask_dict, line[:OrderRef])
                else # If order to be removed IS best ask
                    delete!(Ask_dict, line[:OrderRef]) # Remove order from dictionary
                    # Check if there are more orders in book
                    if isempty(Ask_dict) # If no more orders print NaN
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else # If there ARE more orders in book, update best ask
                        best_ask = findmin(Ask_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, best_ask[1], best_ask[2], NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            else # Order is from Bid
                # Get orderRef of best bid
                ref = findmax(Bid_dict)[2]
                if ref != line[:OrderRef] # If order to be removed is NOT best bid then remove order from dictionary. Don't need to print new best bid
                    delete!(Bid_dict, line[:OrderRef])
                else # If order to be removed IS best bid
                    delete!(Bid_dict, line[:OrderRef]) # Remove order from dictionary
                    # Check if there are more orders in book
                    if isempty(Bid_dict) # If no more orders print NaN
                        temp = (line[:TimeStamp], line[:Date], "BID", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else # If there ARE more orders in book
                        best_bid = findmax(Bid_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "BID", best_bid[1], best_bid[2], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            end
        elseif line[:Type] == "OM" # Check if messeage is an Order Modify
            # Check if the Order Modify is in bid or ask
            if haskey(Ask_dict, line[:OrderRef]) # If true, order is from Ask
                # Get orderRef of best ask
                ref = findmin(Ask_dict)[2]
                # Is the modified order the best ask?
                if ref != line[:OrderRef] # If OM is not best ask
                    # Get the best ask
                    best_ask = findmin(Ask_dict)[1]
                    modified_ask = line[:Price]
                    # Modify order in ask dictionary
                    Ask_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Compare the modified against the best ask
                    if modified_ask < best_ask[1] # Modified is now best ask
                        # Print an event
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, line[:Price], line[:Quantity], NaN, NaN, "")
                        push!(df, temp)
                    end
                else # OM IS best ask
                    # Modify best ask
                    Ask_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Find current best ask after modification
                    new_best_ask = findmin(Ask_dict)[1]
                    # Print an event
                    temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, new_best_ask[1], new_best_ask[2], NaN, NaN, "")
                    push!(df, temp)
                end
            else # order is from Bid
                # Get orderRef of best bid
                ref = findmax(Bid_dict)[2]
                # Is the modified order the best bid?
                if ref != line[:OrderRef] # If OM is not best bid
                    # Get the best bid
                    best_bid = findmax(Bid_dict)[1]
                    modified_bid = line[:Price]
                    # modify order in ask dictionary
                    Bid_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Compare the modified against the best bid
                    if modified_bid > best_bid[1] # Modified is now best bid
                        # Print an event
                        temp = (line[:TimeStamp], line[:Date], "BID", line[:Price], line[:Quantity], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                else # OM IS best bid
                    # Modify best bid
                    Bid_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Find current best bid after modification
                    new_best_bid = findmax(Bid_dict)[1]
                    # Print an event
                    temp = (line[:TimeStamp], line[:Date], "BID", new_best_bid[1], new_best_bid[2], NaN, NaN, NaN, NaN, "")
                    push!(df, temp)
                end
            end
        elseif line[:Type] == "AT" # Check if messeage is a Trade
            # Check if the Order Reference for trade is in bid or ask
            if haskey(Ask_dict, line[:OrderRef]) # If true, order is from Ask; trade is buyer-initiated
                # Print the trade event
                temp = (line[:TimeStamp], line[:Date], "TRADE", NaN, NaN, NaN, NaN, line[:Price], line[:Quantity], "BuyerInitiated")
                push!(df, temp)
                # Minus trade volume from order volume
                vol = Ask_dict[line[:OrderRef]][2] - line[:Quantity]
                # Check if trade finished the order
                if vol != 0 # Trade did not finish entire order
                    # Modify quantity in order
                    Ask_dict[line[:OrderRef]] = (Ask_dict[line[:OrderRef]][1], vol)
                    # Find current best ask after modification
                    new_best_ask = findmin(Ask_dict)[1]
                    # Print updated best ask
                    temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, new_best_ask[1], new_best_ask[2], NaN, NaN, "")
                    push!(df, temp)
                else
                    # Remove the order as it is finished
                    delete!(Ask_dict, line[:OrderRef])
                    # Print updated best ask
                    if isempty(Ask_dict)
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else
                        # Find current best ask after modification
                        new_best_ask = findmin(Ask_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, new_best_ask[1], new_best_ask[2], NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            elseif haskey(Bid_dict, line[:OrderRef]) # Order is from Bid; trade is seller-initiated
                # Print the trade event
                temp = (line[:TimeStamp], line[:Date], "TRADE", NaN, NaN, NaN, NaN, line[:Price], line[:Quantity], "SellerInitiated")
                push!(df, temp)
                # Minus trade volume from order volume
                vol = Bid_dict[line[:OrderRef]][2] - line[:Quantity]
                # Check if trade finished the order
                if vol != 0 # Trade did not finish entire order
                    # Modify quantity in order
                    Bid_dict[line[:OrderRef]] = (Bid_dict[line[:OrderRef]][1], vol)
                    # Find current best bid after modification
                    new_best_bid = findmax(Bid_dict)[1]
                    # Print updated best bid
                    temp = (line[:TimeStamp], line[:Date], "BID", new_best_bid[1], new_best_bid[2], NaN, NaN, NaN, NaN, "")
                    push!(df, temp)
                else
                    # Remove the order as it is finished
                    delete!(Bid_dict, line[:OrderRef])
                    # Print updated best bid
                    if isempty(Bid_dict)
                        temp = (line[:TimeStamp], line[:Date], "BID", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else
                        # Find current best bid after modification
                        new_best_bid = findmax(Bid_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "BID", new_best_bid[1], new_best_bid[2], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            end
        end
    end
    return df
end

## Phase 3 - Extract additional information from order book
function MakeDetailedL1BAT(data::DataFrame) # Function takes in L1 order book and creates the micro-prices and mid-prices, and inter-arrivals
    n = size(data)[1]
    # Initialise micro and mid-price vectors
    micro_price = zeros(n, 1)
    mid_price = zeros(n, 1)
    # Loop through data frame to make micro and mid price vectors
    for i in 1:n
        # Check what event type it is; only want bids and asks
        if data[:EventType][i] == "ASK" # If it is an ask
            # Get current ask and its volume
            current_ask = data[:Ask][i]
            current_ask_vol = data[:AskVol][i]
            # Get index of current bid
            indexof_current_best_bid = findprev(x -> x == "BID", data[:EventType], i)
            if isnothing(indexof_current_best_bid) # No current bids in data; for start of dataset
                # There are no current bids, therefore micro and mid price are NaNs
                micro_price[i] = NaN
                mid_price[i] = NaN
            else
                # Get the index of current bid
                indexof_current_best_bid = indexof_current_best_bid
                # Get current bid and its volume
                current_bid = data[:Bid][indexof_current_best_bid]
                current_bid_vol = data[:BidVol][indexof_current_best_bid]
                # Compute micro and mid price
                micro_price[i] = (current_bid_vol / (current_bid_vol + current_ask_vol)) * current_bid + (current_ask_vol / (current_bid_vol + current_ask_vol)) * current_ask
                mid_price[i] = 0.5*(current_bid + current_ask)
            end
        elseif data[:EventType][i] == "BID" # the event is a Bid
            # Get current bid and its volume
            current_bid = data[:Bid][i]
            current_bid_vol = data[:BidVol][i]
            # Get index of current ask
            indexof_current_best_ask = findprev(x -> x == "ASK", data[:EventType], i)
            if isnothing(indexof_current_best_ask)
                # There are no current asks, therefore micro and mid price are NaNs
                micro_price[i] = NaN
                mid_price[i] = NaN # No current asks in data; for start of dataset
            else
                # Get the index of current ask
                indexof_current_best_ask = indexof_current_best_ask
                # Get current ask and its volume
                current_ask = data[:Ask][indexof_current_best_ask]
                current_ask_vol = data[:AskVol][indexof_current_best_ask]
                # Compute micro and mid price
                micro_price[i] = (current_bid_vol / (current_bid_vol + current_ask_vol)) * current_bid + (current_ask_vol / (current_bid_vol + current_ask_vol)) * current_ask
                mid_price[i] = 0.5*(current_bid + current_ask)
            end
        else
            # The event is a trade. Can't use NaN as that means at least one side of order book is empty
            micro_price[i] = micro_price[i-1]
            mid_price[i] = mid_price[i-1]
        end
    end
    # Initialise vector of inter-arrivals and mid_price change
    τ = fill(NaN, n, 1)
    # Compute inter-arrivals
    TradeInds = findall(x -> x == "TRADE", data[:EventType])
    TradeTimes = data[:TimeStamp][TradeInds]
    inter_arrivals = diff(TradeTimes)
    τ[TradeInds[1:end-1]] = inter_arrivals
    # Initialise DataFrame for detailed L1 Order Book
    new_df = DataFrame(TimeStamp = Float64[], Date = DateTime[], EventType = String[], Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[], Trade = Float64[], TradeVol = Float64[], TradeSign = String[], MicroPrice = Float64[], MidPrice = Float64[], InterArrivals = Float64[])
    for i in 1:n
        temp = (data[i,1], data[i,2], data[i,3], data[i,4], data[i,5], data[i,6], data[i,7], data[i,8], data[i,9], data[i,10], micro_price[i], mid_price[i], τ[i])
        push!(new_df, temp)
    end
    return new_df
end
function GetDetailedL1BAT(ticker, FullData) # Function to streamline the process of phase 2 and 3
    # Read in the ticker
    data_raw = FullData[ticker]
    dates = Date.(data_raw[:,3])
    dates_unique = unique(dates)
    # Initialise dictionary
    ticker_dict = Dict()
    # Loop through each day of the dataset
    @showprogress "Phase 2 & 3 of "*ticker for i in 1:length(dates_unique)
        # Get data for each day
        inds = findall(x -> x == dates_unique[i], dates)
        tempdata = data_raw[inds,:]
        # Make the detailed L1 order book
        tempL1BAT = MakeL1BAT(tempdata)
        tempDetailedL1BAT = MakeDetailedL1BAT(tempL1BAT)
        # Add the day to dictionary
        push!(ticker_dict, Dates.format(dates_unique[i], "yyyy-mm-dd") => tempDetailedL1BAT)
    end
    # Combine into a single master dataframe
    master_df = ticker_dict[Dates.format(dates_unique[1], "yyyy-mm-dd")]
    for i in 2:length(dates_unique)
        master_df = [master_df; ticker_dict[Dates.format(dates_unique[i], "yyyy-mm-dd")]]
    end
    # Write the entire history as a CSV file
    CSV.write("Test Data/A2X/Clean/A2X_Cleaned_"*ticker*".csv", master_df)
end
#---------------------------------------------------------------------------


### 3. Implement cleaning functions
function CleanData() # Function to bring everything together and create
    # Phase 1:
    files = readdir("Test Data/A2X/Raw/") # Find all the files
    filename = "Test Data/A2X/Raw/".*files[2:end]
    Full_data = combineTAQ(filename) # Get data into usable format
    # Phase 2 & 3:
    for i in keys(Full_data) # Loop through each ticker:
        GetDetailedL1BAT(i, Full_data)
    end
end
CleanData()
#---------------------------------------------------------------------------


### 4. Calculate the frequencies of aggressive trades
function AggressiveTradeFrequencies(data)
    buyerData = data[1]; sellerData = data[2] # Seperate into buyer and seller-initiated
    tickers = collect(keys(buyerData)) # Extract the tickers
    aggressiveTradeFrequency = Vector{Float64}() # Initialise vecor of frequencies
    for ticker in tickers
        percent = (sum(isnan.(buyerData[ticker].Impact)) + sum(isnan.(sellerData[ticker].Impact))) / (size(buyerData[ticker], 1) + size(sellerData[ticker], 1)) # Combine buyer and seller-initiated counts
        push!(aggressiveTradeFrequency, percent)
    end
    return aggressiveTradeFrequency
end
aggressiveTradeFrequencies = load("Real Data/A2X/PriceImpact/A2X_PriceImpact.jld")["A2X_PriceImpact"] |> AggressiveTradeFrequencies # Load price impact data
#---------------------------------------------------------------------------
