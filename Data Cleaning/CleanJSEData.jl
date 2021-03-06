### Title: Data Cleaning - JSE
### Authors: Patrick Chang and Ivan Jericevich
### Function: Convert the raw JSE data into a detailed and usable L1LOB format.
### Structure:
# 1. Preliminaries
# 2. Cleaning functions
# 3. Implement cleaning functions
### Strategy:
# Extract data from 2019-01-01 onwards, only pull out data within the continuous trading session, and only keep Automated Trades (AT)
# Compute and add mid-price, micro-price and inter-arrival information
# Classify trades according to the Lee-Ready rule
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, Dates, ProgressMeter
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
#---------------------------------------------------------------------------


### 2. Cleaning Functions
function MakeCleanTAQ(data::DataFrame) # Convert raw data into usable format
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[], EventType = String[], Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[], Trade = Float64[], TradeVol = Float64[])
    # Loop through each day and extract useful data
    @showprogress "Filtering..." for i in 1:length(dates_unique)
        # Obtain single day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Only keep data within continuous trading
        start = DateTime(tempday) + Hour(9)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        tempdata = tempdata[findall(x -> start <= x && x < close, tempdata[:,1]), :]
        # Loop through each line within the day
        for j in 1:size(tempdata)[1]
            # Pull out line item
            line = tempdata[j,:]
            # Check the event type
            if line[:type] == "ASK" # item is an ask
                # Populate with appropriate data
                temp = (line[:times], "ASK", NaN, NaN, line[:value], line[:size], NaN, NaN)
                push!(master_df, temp)
            elseif line[:type] == "BID" # item is a bid
                # Populate with appropriate data
                temp = (line[:times], "BID", line[:value], line[:size], NaN, NaN, NaN, NaN)
                push!(master_df, temp)
            elseif line[:type] == "TRADE" # item is a trade
                # We only want Automated Trades AT
                if line[:condcode] == "AT"
                    # Populate with appropriate data
                    temp = (line[:times], "TRADE", NaN, NaN, NaN, NaN, line[:value], line[:size])
                    push!(master_df, temp)
                end
            end
        end
    end
    return master_df
end
function MakeDetailedTAQ(data::DataFrame) # Get micro-price, mid-price and inter-arrivals information
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    filter!(x-> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[], EventType = String[], Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[], Trade = Float64[], TradeVol = Float64[], MicroPrice = Float64[], MidPrice = Float64[], InterArrivals = Float64[])
    # Loop through each day and extract useful data
    @showprogress "Building..." for j in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[j]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Size of data for the day
        n = size(tempdata)[1]
        # Initialise micro and mid price vectors
        micro_price = zeros(n, 1)
        mid_price = zeros(n, 1)
        # Loop through data frame to make micro and mid price vectors
        for i in 1:n
            # Check what event type it is; only want bids and asks
            if tempdata[:EventType][i] == "ASK"    # if it is an ask
                # Get current ask and its volume
                current_ask = tempdata[:Ask][i]
                current_ask_vol = tempdata[:AskVol][i]
                # Get index of current bid
                indexof_current_best_bid = findlast(x -> x == "BID", tempdata[:EventType][1:i])
                if isnothing(indexof_current_best_bid)    # no current bids in data; for start of dataset
                    # There are no current bids, therefore micro and mid price are NaNs
                    micro_price[i] = NaN
                    mid_price[i] = NaN
                else
                    # Get the index of current bid
                    indexof_current_best_bid = indexof_current_best_bid
                    # Get current bid and its volume
                    current_bid = tempdata[:Bid][indexof_current_best_bid]
                    current_bid_vol = tempdata[:BidVol][indexof_current_best_bid]
                    # Compute micro and mid price
                    micro_price[i] = (current_bid_vol / (current_bid_vol + current_ask_vol)) * current_bid + (current_ask_vol / (current_bid_vol + current_ask_vol)) * current_ask
                    mid_price[i] = 0.5*(current_bid + current_ask)
                end
            elseif tempdata[:EventType][i] == "BID"     # the event is a bid
                # Get current bid and its volume
                current_bid = tempdata[:Bid][i]
                current_bid_vol = tempdata[:BidVol][i]
                # Get index of current ask
                indexof_current_best_ask = findlast(x -> x == "ASK", tempdata[:EventType][1:i])
                if isnothing(indexof_current_best_ask)
                    # There are no current asks, therefore micro and mid price are NaNs
                    micro_price[i] = NaN
                    mid_price[i] = NaN   # No current asks in data; for start of dataset
                else
                    # Get the index of current ask
                    indexof_current_best_ask = indexof_current_best_ask
                    # Get current ask and its volume
                    current_ask = tempdata[:Ask][indexof_current_best_ask]
                    current_ask_vol = tempdata[:AskVol][indexof_current_best_ask]
                    # Compute micro and mid price
                    micro_price[i] = (current_bid_vol / (current_bid_vol + current_ask_vol)) * current_bid + (current_ask_vol / (current_bid_vol + current_ask_vol)) * current_ask
                    mid_price[i] = 0.5*(current_bid + current_ask)
                end
            else
                if i > 2
                    # The event is a trade. Can't use NaN as that means at least one side of order book is empty
                    micro_price[i] = micro_price[i-1]
                    mid_price[i] = mid_price[i-1]
                end
            end
        end
        # Initialise vector of inter-arrivals and mid_price change
        τ = fill(NaN, n, 1)
        # Compute inter-arrivals
        TradeInds = findall(x -> x == "TRADE", tempdata[:EventType])
        TradeTimes = tempdata[:TimeStamp][TradeInds]
        inter_arrivals = diff(datetime2unix.(TradeTimes))
        τ[TradeInds[1:end-1]] = inter_arrivals
        # Push the data into master dataframe
        for i in 1:n
            temp = (tempdata[i,1], tempdata[i,2], tempdata[i,3], tempdata[i,4], tempdata[i,5],
            tempdata[i,6], tempdata[i,7], tempdata[i,8], micro_price[i], mid_price[i], τ[i])
            push!(master_df, temp)
        end
    end
    return master_df
end
function ClassifyTrades(data::DataFrame) # Classify the trades according to the Lee-Ready rule
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    filter!(x-> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[], EventType = String[],
    Bid = Float64[], BidVol = Float64[], Ask = Float64[], AskVol = Float64[],
    Trade = Float64[], TradeVol = Float64[], MicroPrice = Float64[],
    MidPrice = Float64[], InterArrivals = Float64[], TradeSign = String[])
    # Loop through each day and extract useful data
    @showprogress "Classifying..." for k in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[k]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,2])
        # Extract trade values
        TradeValues = tempdata[tradeinds, 7]
        # Initialise inferred classification
        inferredclassification = fill("", size(tempdata)[1], 1)
        # Loop through each trade and classify according to quote rule
        for j in 1:length(tradeinds)
            # Get trade value
            TradeValue = tempdata[tradeinds[j], 7]
            tempdata[1:tradeinds[j], 10]
            # MidQuoteBeforeTrade = tempdata[findlast(!isnan, tempdata[1:tradeinds[j], 10]), 10]
            MidQuoteBeforeTrade = tempdata[tradeinds[j], 10]
            # Perform the logic checks: begin with quote rule
            if TradeValue > MidQuoteBeforeTrade
                # Transaction is higher than mid price => BuyerInitiated
                inferredclassification[tradeinds[j]] = "BuyerInitiated"
            elseif TradeValue < MidQuoteBeforeTrade
                # Transaction is lower than mid price => SellerInitiated
                inferredclassification[tradeinds[j]] = "SellerInitiated"
            elseif TradeValue == MidQuoteBeforeTrade
                # Quote rule failed, go to tick rule
                if j > 1
                    if TradeValues[j] > TradeValues[j-1]
                        # If trade is higher than previous trade => BuyerInitiated
                        inferredclassification[tradeinds[j]] = "BuyerInitiated"
                    elseif TradeValues[j] < TradeValues[j-1]
                        # If trade is lower than previous trade => SellerInitiated
                        inferredclassification[tradeinds[j]] = "SellerInitiated"

                    elseif TradeValues[j] == TradeValues[j-1]
                        # If trades are the same, find last trade that was different
                        indTradeLast = findlast(x -> x != TradeValues[j], TradeValues[1:j])
                        if isnothing(indTradeLast)
                            # No classification, all trades before was the same
                            inferredclassification[tradeinds[j]] = ""
                        else
                            # Compare against last trade that was different
                            TradeLast = TradeValues[indTradeLast]
                            if TradeValues[j] > TradeLast
                                # If trade is higher than previous trade => BuyerInitiated
                                inferredclassification[tradeinds[j]] = "BuyerInitiated"
                            elseif TradeValues[j] < TradeLast
                                # If trade is lower than previous trade => SellerInitiated
                                inferredclassification[tradeinds[j]] = "SellerInitiated"
                            end
                        end
                    end
                else
                    # If first trade can't be classified
                    inferredclassification[tradeinds[j]] = ""
                end
            end
        end
        # Push the data into master dataframe
        for i in 1:size(tempdata)[1]
            temp = (tempdata[i,1], tempdata[i,2], tempdata[i,3], tempdata[i,4], tempdata[i,5],
            tempdata[i,6], tempdata[i,7], tempdata[i,8], tempdata[i,9], tempdata[i,10],
            tempdata[i,11], inferredclassification[i])
            push!(master_df, temp)
        end
    end
    return master_df
end
#---------------------------------------------------------------------------


### 3. Implement cleaning functions
function GetCleanedData(ticker::String)
    data = CSV.read("Test Data/JSE/Raw/JSERAWTAQ"*ticker*".csv") # Read in data
    clean = MakeCleanTAQ(data) # Get data into usable format
    detail = MakeDetailedTAQ(clean) # Make additional information
    classified = ClassifyTrades(detail) # Classify trades using Lee-Ready
    CSV.write("Test Data/JSE/Clean/JSECleanedTAQ"*ticker*".csv", classified) # Write the entire history as a CSV file
end
GetCleanedData("SBK"); GetCleanedData("NPN"); GetCleanedData("SLM"); GetCleanedData("ABG"); GetCleanedData("AGL"); GetCleanedData("BTI"); GetCleanedData("FSR"); GetCleanedData("NED"); GetCleanedData("SHP"); GetCleanedData("SOL")
#---------------------------------------------------------------------------
