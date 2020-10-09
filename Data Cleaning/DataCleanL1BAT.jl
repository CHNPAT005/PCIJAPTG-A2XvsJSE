## Author: Patrick Chang
# Script file to read in the A2X messeages from flat files
# and contructing the L1 Order book with best bid/ask, trades
# and their associated volume

## Phase 2 and 3 of cleaning
#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings

cd("/Users/patrickchang1/HFT2020")

#---------------------------------------------------------------------------
## Phase 2 - Build L1 Order Book
#---------------------------------------------------------------------------
# Function takes in 1 day of cleaned TAQ messages for 1 asset
# and creates the L1 Bid Ask Trade with their associated volumes.
# Trades are also classified into buyer- or seller-initiated.
function MakeL1BAT(data::DataFrame)
    ## Initialise dictionary for all bids and asks
    Bid_dict = Dict()
    Ask_dict = Dict()

    # Initialise DataFrame for L1 Order Book
    df = DataFrame(TimeStamp = Float64[], Date = DateTime[], EventType = String[],
    Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[],
    Trade = Float64[], TradeVol = Float64[], TradeSign = String[])

    # Run through the TAQ messages to create dataframe of L1 BAT
    for i in 1:size(data)[1]
        line = data[i,:]
        if line[:Type] == "OA"  # check if message is an Order Add
            if line[:Side] == "SELL"
                # Check if Dictionary is empty
                if isempty(Ask_dict)    # if empty simply push new ask into dictionary and df
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
                    if new_ask < best_ask   # if new ask is new best, update DataFrame
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, line[:Price], line[:Quantity], NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            else
                # Check if Dictionary is empty
                if isempty(Bid_dict)    # if empty simply push new ask into dictionary and df
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
                    if new_bid > best_bid   # if new bid is new best, update DataFrame
                        temp = (line[:TimeStamp], line[:Date], "BID", line[:Price], line[:Quantity], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            end
        elseif line[:Type] == "OC"  # check if messeage is an Order Cancel
            # Check if the Order Cancel is in bid or ask
            if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask
                # Get orderRef of best ask
                ref = findmin(Ask_dict)[2]
                if ref != line[:OrderRef]
                    # if order to be removed is NOT best ask then remove order from dictionary
                    # don't need to print new best ask
                    delete!(Ask_dict, line[:OrderRef])
                else
                    # if order to be removed IS best ask then remove order from dictionary
                    delete!(Ask_dict, line[:OrderRef])
                    # THEN check if there are more orders in book
                    if isempty(Ask_dict)    # If true print NaN
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else    # If there ARE more orders in book
                        best_ask = findmin(Ask_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, best_ask[1], best_ask[2], NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            else    # Order is from Bid
                # Get orderRef of best bid
                ref = findmax(Bid_dict)[2]
                if ref != line[:OrderRef]
                    # if order to be removed is NOT best bid then remove order from dictionary
                    # don't need to print new best bid
                    delete!(Bid_dict, line[:OrderRef])
                else
                    # if order to be removed IS best bid then remove order from dictionary
                    delete!(Bid_dict, line[:OrderRef])
                    # THEN check if there are more orders in book
                    if isempty(Bid_dict)    # If true print NaN
                        temp = (line[:TimeStamp], line[:Date], "BID", NaN, NaN, NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    else    # If there ARE more orders in book
                        best_bid = findmax(Bid_dict)[1]
                        temp = (line[:TimeStamp], line[:Date], "BID", best_bid[1], best_bid[2], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                end
            end
        elseif line[:Type] == "OM"  # check if messeage is an Order Modify
            # Check if the Order Modify is in bid or ask
            if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask
                # Get orderRef of best ask
                ref = findmin(Ask_dict)[2]
                # Is the modified order the best ask?
                if ref != line[:OrderRef]   # If OM is not best ask
                    # Get the best ask
                    best_ask = findmin(Ask_dict)[1]
                    modified_ask = line[:Price]
                    # modify order in ask dictionary
                    Ask_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Compare the modified against the best ask
                    if modified_ask < best_ask[1]  # Modified is now best ask
                        # Print an event
                        temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, line[:Price], line[:Quantity], NaN, NaN, "")
                        push!(df, temp)
                    end
                else    # OM IS best ask
                    # Modify best ask
                    Ask_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Find current best ask after modification
                    new_best_ask = findmin(Ask_dict)[1]
                    # Print an event
                    temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, new_best_ask[1], new_best_ask[2], NaN, NaN, "")
                    push!(df, temp)
                end
            else    # order is from Bid
                # Get orderRef of best bid
                ref = findmax(Bid_dict)[2]
                # Is the modified order the best bid?
                if ref != line[:OrderRef]   # If OM is not best bid
                    # Get the best bid
                    best_bid = findmax(Bid_dict)[1]
                    modified_bid = line[:Price]
                    # modify order in ask dictionary
                    Bid_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Compare the modified against the best bid
                    if modified_bid > best_bid[1]  # Modified is now best bid
                        # Print an event
                        temp = (line[:TimeStamp], line[:Date], "BID", line[:Price], line[:Quantity], NaN, NaN, NaN, NaN, "")
                        push!(df, temp)
                    end
                else    # OM IS best bid
                    # Modify best bid
                    Bid_dict[line[:OrderRef]] = (line[:Price], line[:Quantity])
                    # Find current best bid after modification
                    new_best_bid = findmax(Bid_dict)[1]
                    # Print an event
                    temp = (line[:TimeStamp], line[:Date], "BID", new_best_bid[1], new_best_bid[2], NaN, NaN, NaN, NaN, "")
                    push!(df, temp)
                end
            end
        elseif line[:Type] == "AT"  # check if messeage is a Trade
            # Check if the Order Reference for trade is in bid or ask
            if haskey(Ask_dict, line[:OrderRef])    # If true, order is from Ask; trade is buyer-initiated
                # Print the trade event
                temp = (line[:TimeStamp], line[:Date], "TRADE", NaN, NaN, NaN, NaN, line[:Price], line[:Quantity], "BuyerInitiated")
                push!(df, temp)
                # Minus trade volume from order volume
                vol = Ask_dict[line[:OrderRef]][2] - line[:Quantity]
                # Check if trade finished the order
                if vol != 0 # trade did not finish entire order
                    # modify quantity in order
                    Ask_dict[line[:OrderRef]] = (Ask_dict[line[:OrderRef]][1], vol)
                    # Find current best ask after modification
                    new_best_ask = findmin(Ask_dict)[1]
                    # Print updated best ask
                    temp = (line[:TimeStamp], line[:Date], "ASK", NaN, NaN, new_best_ask[1], new_best_ask[2], NaN, NaN, "")
                    push!(df, temp)
                else
                    # remove the order as it is finished
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
            elseif haskey(Bid_dict, line[:OrderRef]) # order is from Bid; trade is seller-initiated
                # Print the trade event
                temp = (line[:TimeStamp], line[:Date], "TRADE", NaN, NaN, NaN, NaN, line[:Price], line[:Quantity], "SellerInitiated")
                push!(df, temp)
                # Minus trade volume from order volume
                vol = Bid_dict[line[:OrderRef]][2] - line[:Quantity]
                # Check if trade finished the order
                if vol != 0 # trade did not finish entire order
                    # modify quantity in order
                    Bid_dict[line[:OrderRef]] = (Bid_dict[line[:OrderRef]][1], vol)
                    # Find current best bid after modification
                    new_best_bid = findmax(Bid_dict)[1]
                    # Print updated best bid
                    temp = (line[:TimeStamp], line[:Date], "BID", new_best_bid[1], new_best_bid[2], NaN, NaN, NaN, NaN, "")
                    push!(df, temp)
                else
                    # remove the order as it is finished
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
            # elseif line[:OrderRef] == 0   # trade type 2
            # These are market orders hitting the hidden order book
            # We ignore these for the time being. These are not prevalent
            # so far has only happened for SLM and only 1 of these trades
            end
        end
    end
    return df
end

#---------------------------------------------------------------------------
## Phase 3 - Extracts additional information from Order Book
#---------------------------------------------------------------------------
# Function takes in L1 order book and creates the micro and mid price,
# inter-arrivals and mid price change
function MakeDetailedL1BAT(data::DataFrame)
    n = size(data)[1]
    # Initialise micro and mid price vectors
    micro_price = zeros(n, 1)
    mid_price = zeros(n, 1)
    # Loop through data frame to make micro and mid price vectors
    for i in 1:n
        # Check what event type it is; only want bids and asks
        if data[:EventType][i] == "ASK"    # if it is an ask
            # Get current ask and its volume
            current_ask = data[:Ask][i]
            current_ask_vol = data[:AskVol][i]
            # Get index of current bid
            indexof_current_best_bid = findlast(x -> x == "BID", data[:EventType][1:i])
            if isnothing(indexof_current_best_bid)    # no current bids in data; for start of dataset
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
        elseif data[:EventType][i] == "BID"     # the event is a bid
            # Get current bid and its volume
            current_bid = data[:Bid][i]
            current_bid_vol = data[:BidVol][i]
            # Get index of current ask
            indexof_current_best_ask = findlast(x -> x == "ASK", data[:EventType][1:i])
            if isnothing(indexof_current_best_ask)
                # There are no current asks, therefore micro and mid price are NaNs
                micro_price[i] = NaN
                mid_price[i] = NaN   # no current asks in data; for start of dataset
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
            # The event is a trade
            # Can't use NaN as that means at least one side of order book is empty
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
    new_df = DataFrame(TimeStamp = Float64[], Date = DateTime[], EventType = String[],
    Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[],
    Trade = Float64[], TradeVol = Float64[], TradeSign = String[], MicroPrice = Float64[],
    MidPrice = Float64[], InterArrivals = Float64[])

    for i in 1:n
        temp = (data[i,1], data[i,2], data[i,3], data[i,4], data[i,5], data[i,6],
        data[i,7], data[i,8], data[i,9], data[i,10], micro_price[i],
        mid_price[i], τ[i])
        push!(new_df, temp)
    end
    return new_df
end

#---------------------------------------------------------------------------
# Read in, clean the data and write into dictionary of dataframes

# Only interested in three tickers; these have the most useful trade messages:
# Mr Price (MRP), Standard Bank (SBK), and Sanlam (SLM)

# Function to streamline the process
function GetDetailedL1BAT(ticker)
    # Read in the ticker
    data_raw = CSV.read("A2X/TAQ_full/TAQ_FULL_"*ticker*".csv")
    dates = Date.(data_raw[:,3])
    dates_unique = unique(dates)

    # Initialise dictionary
    ticker_dict = Dict()

    # Loop through each day of the dataset
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Get data for each day
        inds = findall(x -> x == dates_unique[i], dates)
        tempdata = data_raw[inds,:]

        # Make the detailed L1 order book
        tempL1BAT = MakeL1BAT(tempdata)
        tempDetailedL1BAT = MakeDetailedL1BAT(tempL1BAT)

        # Add the day to dictionary
        push!(ticker_dict, Dates.format(dates_unique[i], "yyyy-mm-dd") => tempDetailedL1BAT)

        # Write the day as a CSV file
        CSV.write("A2X/L1BAT/CSV/L1BAT_"*ticker*"_"*Dates.format(dates_unique[i], "yyyy-mm-dd")*".csv", tempDetailedL1BAT)
    end
    save("A2X/L1BAT/JLD/L1BAT_"*ticker*".jld", ticker_dict)
end


# SBK
GetDetailedL1BAT("SBK")
# MRP
GetDetailedL1BAT("MRP")
# SLM
GetDetailedL1BAT("SLM")
# NPN
GetDetailedL1BAT("NPN")

# ARI
GetDetailedL1BAT("ARI")

# AVI
GetDetailedL1BAT("AVI")

# CML
GetDetailedL1BAT("CML")

# HUG
GetDetailedL1BAT("HUG")

# GRT
GetDetailedL1BAT("GRT")

# SNT
GetDetailedL1BAT("SNT")

# APN
GetDetailedL1BAT("APN")
