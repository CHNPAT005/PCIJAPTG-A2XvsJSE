### Title: Construct OHLC Bar Data
### Authors: Patrick Chang and Ivan Jericevich
### Function: Compute the micro-price bar data for NPN along with candlestick plots as well as the transaction bar data combined for all securities in each exchange
### Structure:
# 1. Preliminaries
# 2. Construct transaction bar data combined for all securities to be used for volume seasonality
# 3. Construct micro-price bar data
# 4. Generate candlestick plots
### Interpreting candlestick plot:
# Bar coloured in => open > close
# Bar empty => open < close
# Red => current close < previous close
# Blue => current close > previous close
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, Dates, ProgressMeter, Plots, LaTeXStrings, TimeSeries
#cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()

cd("C:/Users/Ivan/Documents/PCIJAPTG-A2XvsJSE")
NPNJSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQNPN.csv"); NPNA2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_NPN.csv")[:, 2:end]
DataFrames.rename!(NPNA2X, (:Date => :TimeStamp)); select!(NPNJSE, names(NPNA2X))
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
#---------------------------------------------------------------------------


### 2. Construct transaction bar data combined for all securities
function MakeTransactionBars(data::DataFrame, barsize::Integer) # Create the OHLCV + VWAP bar data combined for all securities using transaction data
    # Extract TimeStamp and transaction data
    data = data[:,[1;7;8]]
    # Filter out non-trade times
    data = data[findall(!isnan, data[:,2]),:]
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    #filter!(x -> Date("2019-01-01") <= x <= Date("2019-07-15"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[], Open = Float64[], High = Float64[], Low = Float64[], Close = Float64[], Volume = Float64[], VWAP = Float64[], N = Float64[])
    # Loop through each day and extract useful data
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Only keep data within continuous trading
        start = DateTime(tempday) + Hour(9)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        # Create the bars
        bars = collect(start:Minute(barsize):close)
        # Loop through the bars to create the OHLCV + VWAP bars
        for j in 2:length(bars)
            # Filter out data within the bar
            bardata = tempdata[findall(x -> x > bars[j-1] && x <= bars[j], tempdata[:,1]), :]
            # Check if bar is empty
            if !isempty(bardata)
                # Not empty; make OHLCV + VWAP
                open = bardata[1,2]
                high = maximum(bardata[:,2])
                low = minimum(bardata[:,2])
                close = bardata[end,2]
                volume = sum(bardata[:,3])
                vwap = sum(bardata[:,2] .* bardata[:,3]) / volume
                n = size(bardata)[1]
                # Create data tuple and push into master_df
                temp = (bars[j], open, high, low, close, volume, vwap, n)
                push!(master_df, temp)
            else
                # Empty; push NaNs
                temp = (bars[j], NaN, NaN, NaN, NaN, NaN, NaN, NaN)
                push!(master_df, temp)
            end
        end
    end
    return master_df
end
function CombineTransactionBarData(granularity::Int64, JSE_tickers::Vector{String}, A2X_tickers::Vector{String})
    JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * JSE_tickers[1] * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * A2X_tickers[1] * ".csv")[:, 2:end]
    DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
    barDataJSE = MakeTransactionBars(JSE, granularity); barDataA2X = MakeTransactionBars(A2X, granularity)
    @showprogress "Computing..." for (jseTicker, a2xTicker) in zip(JSE_tickers[2:end], A2X_tickers[2:end])
        JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQ" * jseTicker * ".csv"); A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_" * a2xTicker * ".csv")[:, 2:end]
        DataFrames.rename!(A2X, (:Date => :TimeStamp)); select!(JSE, names(A2X))
        tempJSE = MakeTransactionBars(JSE, granularity); append!(barDataJSE, tempJSE); tempA2X = MakeTransactionBars(A2X, granularity); append!(barDataA2X, tempA2X)
    end
    function WeightedAverage(v, n)
        if sum(.!isnan.(n)) == 0
            return NaN
        else
            indeces = findall(!isnan, n)
            return sum(v[indeces] .* n[indeces]) / sum(n[indeces])
        end
    end
    gdfJSE = groupby(barDataJSE, :TimeStamp); gdfA2X = groupby(barDataA2X, :TimeStamp) # Split
    dfJSE = combine([:Volume, :N] => (v, n) -> (v = WeightedAverage(v, n)), gdfJSE); dfA2X = combine([:Volume, :N] => (v, n) -> (v = WeightedAverage(v, n)), gdfA2X) # Apply and combine
    CSV.write(string("Test Data/JSE/Bar/CombinedTransactionBars", 10, "min.csv"), barDataJSE); CSV.write(string("Test Data/A2X/Bar/CombinedTransactionBars", 10, "min.csv"), barDataA2X)
end
CombineTransactionBarData(10, JSE_tickers, A2X_tickers)
#---------------------------------------------------------------------------


### 3. Construct micro-price bar data
function MakeMicroPriceBars(data::DataFrame, barsize::Integer)
    # Extract TimeStamp and transaction data
    data = data[:,[1;10]]
    # Filter out non-trade times
    data = data[findall(!isnan, data[:,2]),:]
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    #dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[], Open = Float64[], High = Float64[], Low = Float64[], Close = Float64[])
    # Loop through each day and extract useful data
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> Date(x) == tempday && Time(x) > Time(09, 01, 00), data[:,1]), :]
        # Only keep data within continuous trading
        start = DateTime(tempday) + Hour(9)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        # Create the bars
        bars = collect(start:Minute(barsize):close)
        # Loop through the bars to create the OHLCV + VWAP bars
        for j in 2:length(bars)
            # filter out data within the bar
            bardata = tempdata[findall(x-> x>bars[j-1] && x<=bars[j], tempdata[:,1]), :]
            # Check if bar is empty
            if !isempty(bardata)
                # Not empty; make OHLCV + VWAP
                open = bardata[1,2]
                high = maximum(bardata[:,2])
                low = minimum(bardata[:,2])
                close = bardata[end,2]
                # Create data tuple and push into master_df
                temp = (bars[j], open, high, low, close)
                push!(master_df, temp)
            else
                # Empty; push NaNs
                temp = (bars[j], NaN, NaN, NaN, NaN)
                push!(master_df, temp)
            end
        end
    end
    return master_df
end
for interval in [1, 10, 20]
    # JSE bar data
    barDataJSE = MakeMicroPriceBars(NPNJSE, interval)
    CSV.write(string("Test Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"), barDataJSE)
    # A2X bar data
    barDataA2X = MakeMicroPriceBars(NPNA2X, interval)
    CSV.write(string("Test Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"), barDataA2X)
end
#---------------------------------------------------------------------------


### 4. Generate candlestick plots
function candlestick(data::DataFrame, date::String, interval::Int)
    # Get indicies for the date to plot
    inds = findall(x->x==Date(date), Date.(data[:,1]))
    # Filter out data for correct day
    tempdata = data[inds,:]
    tempdata[:, 2:end] = tempdata[:, 2:end] ./ 100
    tempdata[!,1] = Time.(tempdata[:,1])
    # Convert DataFrame into TimeArray object
    ta = TimeArray(tempdata[:,1:5], timestamp = :TimeStamp)
    # Plot the data
    xtick = interval == 1 ? 10 : 1
    p = plot(ta, seriestype = :candlestick, xrotation = 60, color = [:blue,:red], xticks = xtick, foreground_color_grid = :white, xtickfont = 5, dpi = 300)
    xlabel!(p, L"\textrm{2019-07-15}")
    ylabel!(p, L"\textrm{Price [ZAR]}")
    return p
end
for interval in [1, 10, 20]
    # JSE
    NPNJSEMicroBars = CSV.read(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"))
    candlestickJSE = candlestick(NPNJSEMicroBars, "2019-07-12", interval)
    savefig(candlestickJSE, string("Plots/NPNJSEMicroPriceBars", interval, "min.svg"))
    # A2X
    NPNA2XMicroBars = CSV.read(string("Real Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNA2XMicroBars[:,2:5] = NPNA2XMicroBars[:,2:5] ./ 10^5
    candlestickA2X = candlestick(NPNA2XMicroBars, "2019-07-12", interval)
    savefig(candlestickA2X, string("Plots/NPNA2XMicroPriceBars", interval, "min.svg"))
end
#---------------------------------------------------------------------------
