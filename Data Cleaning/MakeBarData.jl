## Author: Patrick Chang
# Script file to read in the cleaned JSE TAQ data and construct
# various types of bar data for candle stick plots
# construct OHLC + VWAP trade bar data & micro- and mid-price
# bar data

#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots
using Statistics, LaTeXStrings, TimeSeries

cd("/Users/patrickchang1/HFT2020")

SBK = CSV.read("JSE/JSECleanedTAQSBK.csv")
NPN = CSV.read("JSE/JSECleanedTAQNPN.csv")

#---------------------------------------------------------------------------
## Functions to make bar data
#---------------------------------------------------------------------------

## Function to create the OHLCV + VWAP bar data using transaction data
# Note that bar size only takes in minutes as the argument
function MakeTransactionBars(data::DataFrame, barsize::Integer)
    # Extract TimeStamp and transaction data
    data = data[:,[1;7;8]]
    # Filter out non-trade times
    data = data[findall(!isnan, data[:,2]),:]

    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    master_df = DataFrame(TimeStamp = DateTime[],
    Open = Float64[], High = Float64[], Low = Float64[], Close = Float64[],
    Volume = Float64[], VWAP = Float64[])
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
            # filter out data within the bar
            bardata = tempdata[findall(x-> x>bars[j-1] && x<=bars[j], tempdata[:,1]),:]
            # Check if bar is empty
            if !isempty(bardata)
                # Not empty; make OHLCV + VWAP
                open = bardata[1,2]
                high = maximum(bardata[:,2])
                low = minimum(bardata[:,2])
                close = bardata[end,2]
                volume = sum(bardata[:,3])
                vwap = sum(bardata[:,2] .* bardata[:,3]) / volume
                # Create data tuple and push into master_df
                temp = (bars[j], open, high, low, close, volume, vwap)
                push!(master_df, temp)
            else
                # Empty; push NaNs
                temp = (bars[j], NaN, NaN, NaN, NaN, NaN, NaN)
                push!(master_df, temp)
            end
        end
    end
    return master_df
end

## Function to create the OHLC bar data using mid-price data
# Note that bar size only takes in minutes as the argument
function MakeMidPriceBars(data::DataFrame, barsize::Integer)
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
        tempdata = data[findall(x -> x == tempday, dates), :]
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

## Function to create the OHLC bar data using mid-price data
# Note that bar size only takes in minutes as the argument
function MakeMicroPriceBars(data::DataFrame, barsize::Integer)
    # Extract TimeStamp and transaction data
    data = data[:,[1;9]]
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
        tempdata = data[findall(x -> x == tempday, dates), :]
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

#---------------------------------------------------------------------------
## Create and write data into flat files
#---------------------------------------------------------------------------

# SBK trades
SBKTradesBars10min = MakeTransactionBars(SBK, 10)
SBKTradesBars1min = MakeTransactionBars(SBK, 1)

# CSV.write("JSE/BarData/SBKTradesBars10min.csv", SBKTradesBars10min)
# CSV.write("JSE/BarData/SBKTradesBars1min.csv", SBKTradesBars1min)

# SBK mid prices
SBKMidPriceBars10min = MakeMidPriceBars(SBK, 10)
SBKMidPriceBars1min = MakeMidPriceBars(SBK, 1)

# CSV.write("JSE/BarData/SBKMidPriceBars10min.csv", SBKMidPriceBars10min)
# CSV.write("JSE/BarData/SBKMidPriceBars1min.csv", SBKMidPriceBars1min)

# SBK micro price
SBKMicroPriceBars10min = MakeMicroPriceBars(SBK, 10)
SBKMicroPriceBars1min = MakeMicroPriceBars(SBK, 1)

# CSV.write("JSE/BarData/SBKMicroPriceBars10min.csv", SBKMicroPriceBars10min)
# CSV.write("JSE/BarData/SBKMicroPriceBars1min.csv", SBKMicroPriceBars1min)

# NPN trades
NPNTradesBars10min = MakeTransactionBars(NPN, 10)
NPNTradesBars1min = MakeTransactionBars(NPN, 1)

# CSV.write("JSE/BarData/NPNTradesBars10min.csv", NPNTradesBars10min)
# CSV.write("JSE/BarData/NPNTradesBars1min.csv", NPNTradesBars1min)

# NPN mid prices
NPNMidPriceBars10min = MakeMidPriceBars(NPN, 10)
NPNMidPriceBars1min = MakeMidPriceBars(NPN, 1)

# CSV.write("JSE/BarData/NPNMidPriceBars10min.csv", NPNMidPriceBars10min)
# CSV.write("JSE/BarData/NPNMidPriceBars1min.csv", NPNMidPriceBars1min)

# NPN micro price
NPNMicroPriceBars10min = MakeMicroPriceBars(NPN, 10)
NPNMicroPriceBars1min = MakeMicroPriceBars(NPN, 1)

# CSV.write("JSE/BarData/NPNMicroPriceBars10min.csv", NPNMicroPriceBars10min)
# CSV.write("JSE/BarData/NPNMicroPriceBars1min.csv", NPNMicroPriceBars1min)

#---------------------------------------------------------------------------
## Plot the data using candlestick plots
#---------------------------------------------------------------------------
# Load the data
SBKTradesBars10min = CSV.read("JSE/BarData/SBKTradesBars10min.csv")
SBKTradesBars1min = CSV.read("JSE/BarData/SBKTradesBars1min.csv")
NPNTradesBars10min = CSV.read("JSE/BarData/NPNTradesBars10min.csv")
NPNTradesBars1min = CSV.read("JSE/BarData/NPNTradesBars1min.csv")

SBKMidPriceBars10min = CSV.read("JSE/BarData/SBKMidPriceBars10min.csv")
SBKMidPriceBars1min = CSV.read("JSE/BarData/SBKMidPriceBars1min.csv")
NPNMidPriceBars10min = CSV.read("JSE/BarData/NPNMidPriceBars10min.csv")
NPNMidPriceBars1min = CSV.read("JSE/BarData/NPNMidPriceBars1min.csv")

SBKMicroPriceBars10min = CSV.read("JSE/BarData/SBKMicroPriceBars10min.csv")
SBKMicroPriceBars1min = CSV.read("JSE/BarData/SBKMicroPriceBars1min.csv")
NPNMicroPriceBars10min = CSV.read("JSE/BarData/NPNMicroPriceBars10min.csv")
NPNMicroPriceBars1min = CSV.read("JSE/BarData/NPNMicroPriceBars1min.csv")

## Function to plot candlestick plot
# Format of date string is "yyyy-mm-dd".
## Interpreting candlestick plot:
# Bar coloured in => open > close
# Bar empty => open < close
# Red => current close < previous close
# Blue => current close > previous close
function candlestick(data::DataFrame, date::String)
    # Get indicies for the date to plot
    inds = findall(x->x==Date(date), Date.(data[:,1]))
    # Filter out data for correct day
    tempdata = data[inds,:]
    tempdata[!,1] = Time.(tempdata[:,1])
    # Convert DataFrame into TimeArray object
    ta = TimeArray(tempdata[:,1:5], timestamp = :TimeStamp)
    # Plot the data
    p = plot(ta, seriestype = :candlestick, dpi = 300)
    f = "\\textrm{"
    xlabel!(p, latexstring("$f$(Date(date))}"))
    ylabel!(p, L"\textrm{Price ZAR [cents]}")
    return p
end
## Function to add VWAP from the bars ontop of the
# candle stick plots. This is specific for transaction data only
function addVWAP(data::DataFrame, date::String; markersize = 3)
    tempdata = data[findall(x->x==Date(date), Date.(data[:,1])),[1;7]]
    plot!(p, Time.(tempdata[:,1]), tempdata[:,2], seriestype = :scatter, markersize = markersize, markershape = :x)
end

## Transaction bars
# SBK 10min
p = candlestick(SBKTradesBars10min, "2019-07-15")
addVWAP(SBKTradesBars10min, "2019-07-15")

# savefig("Assignment2/Plots/SBKTradesBars10min.svg")
# SBK 1min
p = candlestick(SBKTradesBars1min, "2019-07-15")
addVWAP(SBKTradesBars1min, "2019-07-15", markersize = 1)

# savefig("Assignment2/Plots/SBKTradesBars1min.svg")
# NPN 10min
p = candlestick(NPNTradesBars10min, "2019-07-15")
addVWAP(NPNTradesBars10min, "2019-07-15")
# NPN 1min
p = candlestick(NPNTradesBars1min, "2019-07-15")
addVWAP(NPNTradesBars1min, "2019-07-15", markersize = 1)


## MidPrice bars
# SBK 10min
p = candlestick(SBKMidPriceBars10min, "2019-07-15")

# savefig("Assignment2/Plots/SBKMidPriceBars10min.svg")
# SBK 1min
p = candlestick(SBKMidPriceBars1min, "2019-07-15")

# savefig("Assignment2/Plots/SBKMidPriceBars1min.svg")
# NPN 10min
p = candlestick(NPNMidPriceBars10min, "2019-07-15")
# NPN 1min
p = candlestick(NPNMidPriceBars1min, "2019-07-15")

## MicroPrice bars
# SBK 10min
p = candlestick(SBKMicroPriceBars10min, "2019-07-15")

# savefig("Assignment2/Plots/SBKMicroPriceBars10min.svg")
# SBK 1min
p = candlestick(SBKMicroPriceBars1min, "2019-07-15")

# savefig("Assignment2/Plots/SBKMicroPriceBars1min.svg")
# NPN 10min
p = candlestick(NPNMicroPriceBars10min, "2019-07-15")
# NPN 1min
p = candlestick(NPNMicroPriceBars1min, "2019-07-15")
