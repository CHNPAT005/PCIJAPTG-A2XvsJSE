### Title: Extract relevant price-impact information
### Authors: Patrick Chang and Ivan Jericevich
### Function:
### Structure:
# 1. Preliminaries
# 2. Supporting functions for computing direct costs
# 3. Data extraction for price impact
# 4. Implement price-impact extraction functions
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, JLD, Dates, ProgressMeter # , Plots, LaTeXStrings, TimeSeries, Distributions, StatsBase
# cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
cd("C:/Users/Ivan/Documents/PCIJAPTG-A2XvsJSE")
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
#---------------------------------------------------------------------------


### 2. Supporting functions for computing direct costs
# The Ceiling value is in terms of Rands
# JSE values are reported in ZAC
# A2X values are reported in ZAC * 10^5
function JSEDirectCost(TransactionValue; bpsTrade = 0.48, TradeMax = 420.4, bpsSettle = 0.36, SettleMax = 180)
    # Transaction fee
    TradeFee = min(TransactionValue * bpsTrade / 100, TradeMax*10^2)
    SettleFee = min(TransactionValue * bpsSettle / 100, SettleMax*10^2)
    return TradeFee + SettleFee
end
function A2XDirectCost(TransactionValue; bpsTrade = 0.4, TradeMax = 355, bpsSettle = 0.29, SettleMax = 154)
    # Transaction fee
    TradeFee = min(TransactionValue * bpsTrade / 100, TradeMax*10^7)
    SettleFee = min(TransactionValue * bpsSettle / 100, SettleMax*10^7)
    return TradeFee + SettleFee
end
#---------------------------------------------------------------------------


### 3. Data extraction for price impact
function CleanA2X(data::DataFrame) # Further cleaning of the A2X data is necessary before any extraction can occur
    # Extract the dates
    dates = Date.(data[:,2])
    # Get unique dates
    dates_unique = unique(dates)
    # Only keep dates between 2019-01-01 and 2019-07-15
    dates_unique = filter(x-> x >= Date("2019-01-01") && x <= Date("2019-07-15"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(Date = DateTime[], EventType = String[],
    Trade = Float64[], TradeVol = Float64[],
    MidPrice = Float64[], TradeSign = String[])
    # Loop through each day and extract useful data
    for k in 1:length(dates_unique)
        # Extract data from each day
        tempday = dates_unique[k]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,3])
        # Extract trade values
        TradeValues = tempdata[tradeinds, 8]
        # Initialise inferred classification
        inferredclassification = fill("", size(tempdata)[1], 1)
        # Loop through each trade and classify according to quote rule
        for j in 1:length(tradeinds)
            # Get trade value
            TradeValue = tempdata[tradeinds[j], 8]
            tempdata[1:tradeinds[j], 12]
            MidQuoteBeforeTrade = tempdata[findlast(!isnan, tempdata[1:tradeinds[j], 12]), 12]
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
            temp = (tempdata[i,2], tempdata[i,3], tempdata[i,8], tempdata[i,9], tempdata[i,12], inferredclassification[i])
            push!(master_df, temp)
        end
    end
    return master_df
end
function getPriceImpactInfoJSE(data::DataFrame, side::Symbol) # Extract relavent information to create price impact for JSE dataset
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    N = length(dates_unique)
    # Initialise Tⱼ, Volume per day and value traded each day
    Tⱼ = zeros(length(dates_unique), 1)
    Volonday = zeros(length(dates_unique), 1)
    Valueofday = zeros(length(dates_unique), 1)
    # Get Tⱼ for length of dataset
    for i in 1:N
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Get trades for the day
        inds = findall(x->x=="TRADE", tempdata[:,2])
        temptradevols = tempdata[inds,8]
        temptradevalues = tempdata[inds,7]
        # Compute useful info
        Tⱼ[i] = length(inds)
        Volonday[i] = sum(temptradevols)
        Valueofday[i] = sum(temptradevols .* temptradevalues)
    end
    sumTⱼ = sum(Tⱼ)
    AVD = mean(Valueofday)
    # Create master dataframe to store all useful information
    master_df = DataFrame(Impact = Float64[], NormalisedVolume = Float64[], Volume = Float64[], Cost = Float64[], DirectCost = Float64[], Spread = Float64[])
    # Loop through each day and extract useful data
    for i in 1:N
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # If first item is a trade, remove it
        if tempdata[1,2] == "TRADE"
            tempdata = tempdata[2:end,:]
        end
        # Get trades for the day
        inds = findall(x->x=="TRADE", tempdata[:,2])
        # Remove trades with missing classifications
        inds = inds[findall(!ismissing, tempdata[inds,12])]
        # Check if we are getting buy or sell side
        if side == :buy
            signs = tempdata[inds, 12]
            indswewant = findall(x->x=="BuyerInitiated", signs)
        elseif side == :sell
            signs = tempdata[inds, 12]
            indswewant = findall(x->x=="SellerInitiated", signs)
        else
            return println("Pick correct option for side with either :buy or :sell")
        end
        # Get indicies that we want
        inds = inds[indswewant]
        # Loop through the trades for each day
        for j in 1:length(inds)
            # Compute the relavent information
            MPbefore = tempdata[inds[j]-1, 10]; MPafter = tempdata[inds[j]+1, 10]
            Δp = log(MPafter) - log(MPbefore)
            # Δp = log(tempdata[inds[j]+1, 10]) - log(tempdata[inds[j]-1, 10])
            vol = tempdata[inds[j], 8]
            normvol = (vol / Volonday[i]) * (sumTⱼ / N)
            TradePrice = tempdata[inds[j], 7]
            DirectCost = JSEDirectCost(TradePrice * vol) / vol
            s = abs(TradePrice - MPbefore)
            # Push to master_df; add correct sign to price change
            if side == :buy
                Δc = log(MPafter + DirectCost + s) - log(MPbefore)
                Δd = log(MPbefore + DirectCost) - log(MPbefore)
                Δs = log(MPbefore + s) - log(MPbefore)
                temp = (Δp, normvol, vol, Δc, Δd, Δs)
                push!(master_df, temp)
            elseif side == :sell
                Δc = log(MPafter - DirectCost - s) - log(MPbefore)
                Δd = log(MPbefore + DirectCost) - log(MPbefore)
                Δs = log(MPbefore + s) - log(MPbefore)
                temp = (-Δp, normvol, vol, -Δc, Δd, Δs)
                push!(master_df, temp)
            end
        end
    end
    return master_df, AVD
end
function getPriceImpactInfoA2X(data::DataFrame, side::Symbol) # Extract relavent information to create price impact for A2X dataset
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    N = length(dates_unique)
    # Initialise Tⱼ, Volume per day and value traded each day
    Tⱼ = zeros(length(dates_unique), 1)
    Volonday = zeros(length(dates_unique), 1)
    Valueofday = zeros(length(dates_unique), 1)
    # Get Tⱼ for length of dataset
    for i in 1:N
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Get trades for the day
        inds = findall(x->x=="TRADE", tempdata[:,2])
        temptradevols = tempdata[inds,4]
        temptradevalues = tempdata[inds,3]
        # Compute useful info
        Tⱼ[i] = length(inds)
        Volonday[i] = sum(temptradevols)
        Valueofday[i] = sum(temptradevols .* temptradevalues)
    end
    sumTⱼ = sum(Tⱼ)
    ADV = mean(Valueofday)
    # Create master dataframe to store all useful information
    master_df = DataFrame(Impact = Float64[], NormalisedVolume = Float64[], Volume = Float64[], Cost = Float64[], DirectCost = Float64[], Spread = Float64[])
    # Loop through each day and extract useful data
    for i in 1:N
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # tempdata = tempdata[findall(!isnan, tempdata[:,5]),:]
        # Get trades for the day
        inds = findall(x->x=="TRADE", tempdata[:,2])
        # Check if we are getting buy or sell side
        if side == :buy
            signs = tempdata[inds, 6]
            indswewant = findall(x->x=="BuyerInitiated", signs)
        elseif side == :sell
            signs = tempdata[inds, 6]
            indswewant = findall(x->x=="SellerInitiated", signs)
        else
            return println("Pick correct option for side with either :buy or :sell")
        end
        # Get indicies that we want
        inds = inds[indswewant]
        # Loop through the trades for each day
        for j in 1:length(inds)
            # Compute the relavent information
            MPbefore = tempdata[inds[j]-1, 5]; MPafter = tempdata[inds[j]+1, 5]
            Δp = log(MPafter) - log(MPbefore)
            # Δp = log(tempdata[inds[j]+1, 5]) - log(tempdata[inds[j]-1, 5])
            vol = tempdata[inds[j], 4]
            normvol = (vol / Volonday[i]) * (sumTⱼ / N)
            TradePrice = tempdata[inds[j], 3]
            DirectCost = A2XDirectCost(TradePrice * vol) / vol
            s = abs(TradePrice - MPbefore)
            # Push to master_df; add correct sign to price change
            if side == :buy
                Δc = log(MPafter + DirectCost + s) - log(MPbefore)
                Δd = log(MPbefore + DirectCost) - log(MPbefore)
                Δs = log(MPbefore + s) - log(MPbefore)
                temp = (Δp, normvol, vol, Δc, Δd, Δs)
                push!(master_df, temp)
            elseif side == :sell
                Δc = log(MPafter - DirectCost - s) - log(MPbefore)
                Δd = log(MPbefore + DirectCost) - log(MPbefore)
                Δs = log(MPbefore + s) - log(MPbefore)
                temp = (-Δp, normvol, vol, -Δc, Δd, Δs)
                push!(master_df, temp)
            end
        end
    end
    return master_df, ADV
end
function PriceImpactInfoA2X(ticker::Vector) # Read in the A2X data, create seller- and buyer- initiated dictionaries containing all the relavent price impact information for the various tickers
    # Create the relavent dictionaries
    BuyerDict = Dict(); SellerDict = Dict(); ADVDict = Dict()
    # Read in the data
    @showprogress "Computing..." for i in 1:length(ticker)
        # Read in the data
        data = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_"*ticker[i]*".csv")
        # Clean the A2X data and return dataframe of relavent information
        CleanedData = CleanA2X(data)
        # Extract price impact information
        Buyer = getPriceImpactInfoA2X(CleanedData, :buy)
        ADV = Buyer[2]
        BuyerInitiated = Buyer[1]
        SellerInitiated = getPriceImpactInfoA2X(CleanedData, :sell)[1]
        # Push the price impact information into their appropriate dictionaries
        push!(BuyerDict, ticker[i] => BuyerInitiated)
        push!(SellerDict, ticker[i] => SellerInitiated)
        push!(ADVDict, ticker[i] => ADV)
    end
    return BuyerDict, SellerDict, ADVDict
end
function PriceImpactInfoJSE(ticker::Vector) # Read in the JSE data, create seller- and buyer- initiated dictionaries containing all the relavent price impact information for the various tickers
    # Create the relavent dictionaries
    BuyerDict = Dict(); SellerDict = Dict(); ADVDict = Dict()
    # Read in the data
    @showprogress "Computing..." for i in 1:length(ticker)
        # Read in the data
        data = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ"*ticker[i]*".csv")
        # Extract price impact information
        Buyer = getPriceImpactInfoJSE(data, :buy)
        ADV = Buyer[2]
        BuyerInitiated = Buyer[1]
        SellerInitiated = getPriceImpactInfoJSE(data, :sell)[1]
        # Push the price impact information into their appropriate dictionaries
        push!(BuyerDict, ticker[i] => BuyerInitiated)
        push!(SellerDict, ticker[i] => SellerInitiated)
        push!(ADVDict, ticker[i] => ADV)
    end
    return BuyerDict, SellerDict, ADVDict
end
#---------------------------------------------------------------------------


# 4. Implement price-impact extraction functions
A2X_PriceImpact = PriceImpactInfoA2X(A2X_tickers); save("Test Data/A2X/Price Impact/A2X_PriceImpact.jld", "A2X_PriceImpact", A2X_PriceImpact)
JSE_PriceImpact = PriceImpactInfoJSE(JSE_tickers); save("Test Data/JSE/Price Impact/JSE_PriceImpact.jld", "JSE_PriceImpact", JSE_PriceImpact)
#---------------------------------------------------------------------------
