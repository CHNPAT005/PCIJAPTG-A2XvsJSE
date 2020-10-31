### Title: Seasonality
### Authors: Patrick Chang and Ivan Jericevich
### Function: Investigate the intraday order-book seasonality of JSE and A2X
### Structure:
# 1. Preliminaries
# 2. Aggregate volume across the day normalised by daily volume
# 3. Average absolute intraday returns normalised by absolute intraday returns
# 4. Average spread normalised by the average daily spread
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
#---------------------------------------------------------------------------


### 2. Aggregate volume across the day normalised by daily volume
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
JSETradesBars10min = CSV.File("Test Data/JSE/Bar/CombinedTransactionBars10min.csv") |> DataFrame!; A2XTradesBars10min = CSV.File("Test Data/A2X/Bar/CombinedTransactionBars10min.csv") |> DataFrame!
# Compute the volume curve and save data
JSEAveNormVol = AveNormVol(JSETradesBars10min); A2XAveNormVol = AveNormVol(A2XTradesBars10min); save("Computed Data/VolumeSeasonality.jld", "JSEVolumeSeasonality", JSEAveNormVol, "A2XVolumeSeasonality", A2XAveNormVol)
# Load the pre-computed data
volumeSeasonality = load("Computed Data/VolumeSeasonality.jld"); JSEVolumeSeasonality = volumeSeasonality["JSEVolumeSeasonality"]; A2XVolumeSeasonality = volumeSeasonality["A2XVolumeSeasonality"]
plot(JSEVolumeSeasonality[1], JSEVolumeSeasonality[2], seriestype = :bar, label = "", fillcolor = :red, dpi = 300, legend = :topleft, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Volume}"); savefig("Figures/JSEVolumeSeasonality.pdf")
plot(A2XVolumeSeasonality[1], A2XVolumeSeasonality[2], seriestype = :bar, label = "", fillcolor = :blue, dpi = 300, legend = :topleft, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Volume}"); savefig("Figures/A2XVolumeSeasonality.svg")
#---------------------------------------------------------------------------


### 3. Average absolute intraday returns normalised by absolute intraday returns
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
        tempdata = tempdata[findall(!isnan, tempdata[:, :Trade]), :]
        # Pull out the item to make returns with
        prices = tempdata[:, :Trade]
        ret = diff(log.(prices))
        # Append to vector of returns
        append!(Master_df, DataFrame(TimeStamp = tempdata[2:end, 1], Return = ret))
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
JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * JSE_tickers[1] * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * A2X_tickers[1] * ".csv")[:, 2:end]
DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
returnJSE = getReturns(JSE); returnA2X = getReturns(A2X)
for (jseTicker, a2xTicker) in zip(JSE_tickers[2:end], A2X_tickers[2:end])
    JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * jseTicker * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * a2xTicker * ".csv")[:, 2:end]
    # Ensure that the column names and ordering are the same
    DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
    tempJSE = getReturns(JSE); tempA2X = getReturns(A2X)
    # Concactenate tickers
    append!(returnJSE, tempJSE); append!(returnA2X, tempA2X)
end
# Compute the return curve and save data
JSEReturnSeasonality = AveAbsRet(returnJSE, 10); A2XReturnSeasonality = AveAbsRet(returnA2X, 10)
save("Computed Data/ReturnSeasonality.jld", "JSEReturnSeasonality", JSEReturnSeasonality, "A2XReturnSeasonality", A2XReturnSeasonality)
# Load the pre-computed data
returnSeasonality = load("Computed Data/ReturnSeasonality.jld"); JSEReturnSeasonality = returnSeasonality["JSEReturnSeasonality"]; A2XReturnSeasonality = returnSeasonality["A2XReturnSeasonality"]
plot(JSEReturnSeasonality[1], JSEReturnSeasonality[2], seriestype = :bar, label = "", fillcolor = :red, dpi = 300, legend = :topright, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Absolute Returns}"); savefig("Figures/JSEReturnSeasonality.pdf")
plot(A2XReturnSeasonality[1], A2XReturnSeasonality[2], seriestype = :bar, label = "", fillcolor = :blue, dpi = 300, legend = :topright, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Absolute Returns}"); savefig("Figures/A2XReturnSeasonality.pdf")
#---------------------------------------------------------------------------


### 4. Average spread normalised by the average daily spread
function getSpread(data::DataFrame)
    # Initialise Master dataframe to store the timestamp + spread
    Master_df = DataFrame(TimeStamp = DateTime[], Spread = Float64[])
    # Get unique dates and start from 2019-01-01 onwards
    filter!(x-> Date("2019-01-01") <= Date(x.TimeStamp) <= Date("2019-07-15"), data)
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
JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * JSE_tickers[1] * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * A2X_tickers[1] * ".csv")[:, 2:end]
DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
spreadJSE = getSpread(JSE); spreadA2X = getSpread(A2X)
for (jseTicker, a2xTicker) in zip(JSE_tickers[2:end], A2X_tickers[2:end])
    JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * jseTicker * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * a2xTicker * ".csv")[:, 2:end]
    # Ensure that the column names and ordering are the same
    DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
    tempJSE = getSpread(JSE); tempA2X = getSpread(A2X)
    # Concactenate tickers
    append!(spreadJSE, tempJSE); append!(spreadA2X, tempA2X)
end
# Compute the return curve and save data
JSESpreadSeasonality = AveSpread(spreadJSE, 10); A2XSpreadSeasonality = AveSpread(spreadA2X, 10)
save("Computed Data/SpreadSeasonality.jld", "JSESpreadSeasonality", JSESpreadSeasonality, "A2XSpreadSeasonality", A2XSpreadSeasonality)
# Load the pre-computed data
spreadSeasonality = load("Computed Data/SpreadSeasonality.jld"); JSESpreadSeasonality = spreadSeasonality["JSESpreadSeasonality"]; A2XSpreadSeasonality = spreadSeasonality["A2XSpreadSeasonality"]
plot(JSESpreadSeasonality[1], JSESpreadSeasonality[2], seriestype = :bar, label = "", dpi = 300, legend = :topright, fillcolor = :red, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Spread}"); savefig("Figures/JSESpreadSeasonality.svg")
plot(A2XSpreadSeasonality[1], A2XSpreadSeasonality[2], seriestype = :bar, label = "", fillcolor = :blue, dpi = 300, legend = :topright, xlabel = L"\textrm{Time of day}", ylabel = L"\textrm{Normalised Spread}"); savefig("Figures/A2XSpreadSeasonality.svg")
#---------------------------------------------------------------------------
