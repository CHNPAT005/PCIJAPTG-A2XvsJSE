### Title: Bar Data
### Authors: Patrick Chang and Ivan Jericevich
### Function: Generate micro-price bar data and candlestick plots
### Structure:
# 1. Preliminaries
# 2. Construct micro-price bar data
# 3. Generate candlestick plots
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings, TimeSeries
cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")
NPNJSE = CSV.read("Real Data/JSE/Cleaned/JSECleanedTAQNPN.csv")
NPNA2X = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")[:, 2:end]
DataFrames.rename!(NPNA2X, (:Date => :TimeStamp))
select!(NPNJSE, names(NPNA2X))
#---------------------------------------------------------------------------


### 2. Construct micro-price bar data
function MakeMicroPriceBars(data::DataFrame, barsize::Integer)
    # Extract TimeStamp and transaction data
    data = data[:,[1;10]]
    # Filter out non-trade times
    data = data[findall(!isnan, data[:,2]),:]
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[],
    Open = Float64[], High = Float64[], Low = Float64[], Close = Float64[])
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
            bardata = tempdata[findall(x-> x>bars[j-1] && x<=bars[j], tempdata[:,1]),:]
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
for interval in [1, 10, 20, 30]
    # JSE bar data
    barDataJSE = MakeMicroPriceBars(NPNJSE, interval)
    CSV.write(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"), barDataJSE)
    # A2X bar data
    barDataA2X = MakeMicroPriceBars(NPNA2X, interval)
    CSV.write(string("Real Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"), barDataA2X)
end
#---------------------------------------------------------------------------


### 3. Generate candlestick plots
function candlestick(data::DataFrame, date::String, interval::Int)
    ## Interpreting candlestick plot:
    # Bar coloured in => open > close
    # Bar empty => open < close
    # Red => current close < previous close
    # Blue => current close > previous close
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
    p = plot(ta, seriestype = :candlestick, xrotation = 60, color = [:blue,:red], xticks = xtick, foreground_color_grid = :white)
    f = "\\textrm{"
    xlabel!(p, L"\\textrm{2019-07-15}")
    ylabel!(p, L"\\textrm{Price ZAR [cents]}")
    return p
end
for interval in [1, 10, 20, 30]
    # JSE
    NPNJSEMicroBars = CSV.read(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"))
    candlestickJSE = candlestick(NPNJSEMicroBars, "2019-07-15", interval)
    savefig(candlestickJSE, string("Plots/NPNJSEMicroPriceBars", interval, "min.pdf"))
    # A2X
    NPNA2XMicroBars = CSV.read(string("Real Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"))
    candlestickA2X = candlestick(NPNA2XMicroBars, "2019-07-15", interval)
    savefig(candlestickA2X, string("Plots/NPNA2XMicroPriceBars", interval, "min.pdf"))
end
#---------------------------------------------------------------------------
