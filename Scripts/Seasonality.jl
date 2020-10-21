## Author: Patrick Chang
# Script file to investigate the order-book seasonality of JSE
# Here we split the data into two sections based on daylight savings
# which occured on 31 March 2019
# We are looking at:
#   1: Aggregate volume across the day normalised by daily volume
#   2: Average absolute intraday returns normalised by absolute intraday returns
#   3: Average spread normalised by the average daily spread
#---------------------------------------------------------------------------
## Preamble
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots
using Statistics, LaTeXStrings, TimeSeries
cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]
A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]

function CombineTransactionBarData(exchange::String, granularity::Int64, tickers::Vector{String})
    TradesBars = CSV.read(string("Real Data/", exchange, "/Bar/", tickers[1], "TradesBars", granularity, "min.csv"))# |> y -> filter(x -> !isnan(x.N), y)
    @showprogress "Computing..." for ticker in tickers[2:end]
        tempBars = CSV.read(string("Real Data/", exchange, "/Bar/", ticker, "TradesBars", granularity, "min.csv"))# |> y -> filter(x -> !isnan(x.N), y)
        append!(TradesBars, tempBars)
    end
    function WeightedAverage(v, n)
        if sum(.!isnan.(n)) == 0
            return NaN
        else
            indeces = findall(.!isnan.(n))
            return sum(v[indeces] .* n[indeces]) / sum(n[indeces])
        end
    end
    #sort!(TradesBars10minBDT, ); sort!(TradesBars10minBDST)
    gdf = groupby(TradesBars, :TimeStamp)
    df = combine([:Volume, :N] => (v, n) -> (v = WeightedAverage(v, n)), gdf) # Apply and combine
    return df
end
JSETradesBars1min = CombineTransactionBarData("JSE", 1, JSE_tickers); JSETradesBars10min = CombineTransactionBarData("JSE", 10, JSE_tickers)
A2XTradesBars1min = CombineTransactionBarData("A2X", 1, A2X_tickers); A2XTradesBars10min = CombineTransactionBarData("A2X", 10, A2X_tickers)
#---------------------------------------------------------------------------
# 1: Aggregate volume across the day normalised by daily volume
#---------------------------------------------------------------------------
## Function to compute the average volume across the day normalised by daily volume
function AveNormVol(data::DataFrame)
    # Extract the dates of the data
    dates = Date.(data[:,1])
    dates_unique = unique(dates)
    # Number of bins of bars each day
    nbins = length(findall(x->x==dates_unique[1], dates))
    # Create a matrix to store the normalised volume of each bin over the various days
    NormData = zeros(nbins, length(dates_unique))
    # Loop through each day of data
    for i in 1:length(dates_unique)
        # Extract the data for the day
        tempdata = data[findall(x->x==dates_unique[i], dates),:]
        # Get total vol for the day
        DayVol = sum(tempdata[.!isnan.(tempdata[:, 2]),2])
        # Loop through each bin of the day
        for j in 1:nbins
            NormData[j, i] = tempdata[j, 2] / DayVol
        end
    end
    # Initialise vector of aggregated normalised volume
    AggregatedNormData = zeros(nbins, 1)
    # Compute the average norm vol for each bin and removing any empty bins
    for i in 1:nbins
        AggregatedNormData[i] = mean(filter(!isnan, NormData[i,:]))
    end
    times = Time.(data[1:nbins,1])
    return times, AggregatedNormData
end

JSEAveNormVol = AveNormVol(JSETradesBars10min); A2XAveNormVol = AveNormVol(A2XTradesBars10min)
plot(A2XAveNormVol[1], A2XAveNormVol[2], seriestype = :bar, label = L"\textrm{A2X}", fillcolor = :blue, dpi = 300, legend = :topleft)
plot!(JSEAveNormVol[1], JSEAveNormVol[2], seriestype = :bar, label = L"\textrm{JSE}", fillcolor = :red)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Volume}")
# savefig("Assignment2/Plots/SBK_AveNormVol.svg")



#---------------------------------------------------------------------------
# 2: Average absolute intraday returns normalised by absolute intraday returns
#---------------------------------------------------------------------------
## Function to compute the average absolute intraday returns normalised by
# the average daily absolute intraday returns. We use VWAP prices from each bar.
function getReturns(data::DataFrame)
    # Create master dataframe to store all useful information
    Master_df = DataFrame(TimeStamp = DateTime[], Return = Float64[])
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    filter!(x-> Date("2019-01-01") <= x <= Date("2019-07-15"), dates_unique)
    # Loop through each day and piece the returns from there
    @showprogress "Computing..." for k in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[k]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Pull out the item to make returns with
        prices = tempdata[:, :MicroPrice]
        # Returns are price fluctuations
        indeces = findall(!isnan, prices)
        ret = diff(log.(prices[indeces]))
        # Append to vector of returns
        append!(Master_df, DataFrame(TimeStamp = tempdata[indeces[2:end], 1], Return = ret))
    end
    return Master_df
end

function AveAbsRet(data::DataFrame, barsize::Integer)
    # Extract the dates of the data
    dates = Date.(data[:,1])
    dates_unique = unique(dates)
    # Number of bins of bars each day
    nbins = length(collect(Time(9,0):Minute(barsize):Time(16,50))) - 1
    # Create a matrix to store the normalised spread of each bin over the various days
    NormData = zeros(nbins, length(dates_unique))
    # Loop through each day of data
    for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x->x == tempday, dates),:]
        # Only keep data within continuous trading and remove the first minute corresponding to the auction call
        start = DateTime(tempday) + Hour(9) + Minute(10)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        # Create the bars
        bars = vcat(DateTime(tempday) + Hour(9) + Minute(1), collect(start:Minute(barsize):close))
        # Get average absolute intraday bar returns
        AveAbsDailyRet = mean(abs.(tempdata[:,2]))
        # Loop through each bin of the day
        for j in 1:(length(bars)-1)
            # filter out data within the bar
            bardata = tempdata[findall(x-> x>bars[j] && x<=bars[j+1], tempdata[:,1]),:]
            # Average spread in the bar normalised by average spread for the day
            NormData[j, i] = mean(abs.(bardata[:,2])) / AveAbsDailyRet
        end
    end
    # Initialise vector of aggregated absolute intraday returns
    AggregatedNormData = zeros(nbins, 1)
    # Compute the average absolute intraday returns for each bin and removing any empty bins
    for i in 1:nbins
        AggregatedNormData[i] = mean(filter(!isnan, NormData[i,:]))
    end
    times = collect(Time(9,10):Minute(barsize):Time(16,50))
    return times, AggregatedNormData
end
# Read full cleaned data files for each ticker and concactenate spreads calculated tick-by-tick
JSE = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/JSECleanedTAQ" * JSE_tickers[1] * ".csv")
A2X = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/A2X_Cleaned_" * A2X_tickers[1] * ".csv")[:, 2:end]
DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
returnJSE = getReturns(JSE); returnA2X = getReturns(A2X)
for (jseTicker, a2xTicker) in zip(JSE_tickers[2:end], A2X_tickers[2:end])
    JSE = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/JSECleanedTAQ" * jseTicker * ".csv")
    A2X = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/A2X_Cleaned_" * a2xTicker * ".csv")[:, 2:end]
    # Ensure that the column names and ordering are the same
    DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
    tempJSE = getReturns(JSE); tempA2X = getReturns(A2X)
    # Concactenate tickers
    append!(returnJSE, tempJSE); append!(returnA2X, tempA2X)
end
# save("returnJSE.jld", "returnJSE", returnJSE); save("returnA2X.jld", "returnA2X", returnA2X)# returnJSE = load("returnJSE.jld")["returnJSE"]; returnA2X = load("returnA2X.jld")["returnA2X"]
JSEReturnSeasonality = AveAbsRet(returnJSE, 10); A2XReturnSeasonality = AveAbsRet(returnA2X, 10)
save("ReturnSeasonality.jld", "JSEReturnSeasonality", JSEReturnSeasonality, "A2XReturnSeasonality", A2XReturnSeasonality)
returnSeasonality = load("Computed Data/ReturnSeasonality.jld"); JSEReturnSeasonality = returnSeasonality["JSEReturnSeasonality"]; A2XReturnSeasonality = returnSeasonality["A2XReturnSeasonality"]
# JSEBuy
plot(JSEReturnSeasonality[1], JSEReturnSeasonality[2], seriestype = :bar, label = L"\textrm{BDST}", fillcolor = :red, dpi = 300, legend = :topright)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Absolute Returns}")
# A2X
plot(A2XReturnSeasonality[1], A2XReturnSeasonality[2], seriestype = :bar, label = L"\textrm{BDST}", fillcolor = :blue, dpi = 300, legend = :topright)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Absolute Returns}")

#---------------------------------------------------------------------------
# 3: Average spread normalised by the average daily spread
#---------------------------------------------------------------------------
## Function to get the spread data
function getSpread(data::DataFrame)
    # Initialise Master dataframe to store the timestamp + spread
    Master_df = DataFrame(TimeStamp = DateTime[], Spread = Float64[])
    # Get unique dates and start from 2019-01-01 onwards
    filter!(x-> Date("2019-01-01") <= Date(x) <= Date("2019-07-15"), data[:, 1])
    # Loop through each item of data and push data into Master_df
    @showprogress "Computing..." for i in 1:size(data)[1]
        line = data[i,:]
        # If there is no mid-price then we don't compute spread
        if !isnan(line[11])
            # Check if current update is best bid or ask
            if line[2] == "BID" # If Bid
                # Get Spread = 2* (midprice - best bid)
                spread = 2 * abs(line[11] - line[3])
                # Push data into Master_df
                push!(Master_df, (line[1], spread))
            elseif line[2] == "ASK" # If ask
                # Get spread: 2 * (best ask - midprice)
                spread = 2 * abs(line[5] - line[11])
                # Push data into Master_df
                push!(Master_df, (line[1], spread))
            end # We don't consider trades when making the spread
        end
    end
    return Master_df
end
## Function to compute the average spread normalised by the average daily spread
function AveSpread(data::DataFrame, barsize::Integer)
    # Extract the dates of the data
    dates = Date.(data[:,1])
    dates_unique = unique(dates)
    # Number of bins of bars each day
    nbins = length(collect(Time(9,0):Minute(barsize):Time(16,50))) - 1
    # Create a matrix to store the normalised spread of each bin over the various days
    NormData = zeros(nbins, length(dates_unique))
    # Loop through each day of data
    for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x->x == tempday, dates),:]
        # Only keep data within continuous trading and remove the first minute corresponding to the auction call
        start = DateTime(tempday) + Hour(9) + Minute(10)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        # Create the bars
        bars = vcat(DateTime(tempday) + Hour(9) + Minute(1), collect(start:Minute(barsize):close))
        # Get average spread for the day
        AveDailySpread = mean(tempdata[:,2])
        # Loop through each bin of the day
        for j in 1:(length(bars)-1)
            # filter out data within the bar
            bardata = tempdata[findall(x-> x>bars[j] && x<=bars[j+1], tempdata[:,1]),:]
            # Average spread in the bar normalised by average spread for the day
            NormData[j, i] = mean(bardata[:,2]) / AveDailySpread
        end
    end
    # Initialise vector of aggregated normalised spread
    AggregatedNormData = zeros(nbins, 1)
    # Compute the average norm spread for each bin and removing any empty bins
    for i in 1:nbins
        AggregatedNormData[i] = mean(filter(!isnan, NormData[i,:]))
    end
    times = collect(Time(9,10):Minute(barsize):Time(16,50))
    return times, AggregatedNormData
end
# Read full cleaned data files for each ticker and concactenate spreads calculated tick-by-tick
JSE = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/JSECleanedTAQ" * JSE_tickers[1] * ".csv")
A2X = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/A2X_Cleaned_" * A2X_tickers[1] * ".csv")[:, 2:end]
DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
spreadJSE = getSpread(JSE); spreadA2X = getSpread(A2X)
for (jseTicker, a2xTicker) in zip(JSE_tickers[2:end], A2X_tickers[2:end])
    JSE = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/JSECleanedTAQ" * jseTicker * ".csv")
    A2X = CSV.read("C:/Users/Ivan/University of Cape Town/Patrick Chang - JSE Raw Data/A2X_Cleaned_" * a2xTicker * ".csv")[:, 2:end]
    # Ensure that the column names and ordering are the same
    DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
    tempJSE = getSpread(JSE); tempA2X = getSpread(A2X)
    # Concactenate tickers
    append!(spreadJSE, tempJSE); append!(spreadA2X, tempA2X)
end
# spreadJSE = load("spreadJSE.jld")["spreadJSE"]; spreadA2X = load("spreadA2X.jld")["spreadA2X"] # save("spreadJSE.jld", "spreadJSE", spreadJSE); save("spreadA2X.jld", "spreadA2X", spreadA2X)
JSESpreadSeasonality = AveSpread(spreadJSE, 10); A2XSpreadSeasonality = AveSpread(spreadA2X, 10)
spreadSeasonality = load("Computed Data/SpreadSeasonality.jld"); JSESpreadSeasonality = spreadSeasonality["JSESpreadSeasonality"]; A2XSpreadSeasonality = spreadSeasonality["A2XSpreadSeasonality"] # save("SpreadSeasonality.jld", "JSESpreadSeasonality", JSESpreadSeasonality, "A2XSpreadSeasonality", A2XSpreadSeasonality)
# JSE
plot(JSESpreadSeasonality[1], JSESpreadSeasonality[2], seriestype = :bar, label = L"\textrm{JSE}", dpi = 300, legend = :topright, fillcolor = :red)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Spread}")
# A2X
plot(A2XSpreadSeasonality[1], A2XSpreadSeasonality[2], seriestype = :bar, label = L"\textrm{A2X}", fillcolor = :blue, dpi = 300, legend = :topright)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Spread}")
# savefig("Assignment2/Plots/SBK_Spread.svg")
