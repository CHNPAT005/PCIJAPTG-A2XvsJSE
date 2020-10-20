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
cd("C:/Users/Ivan/Documents/PCIJAPTG-A2XvsJSE")
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]
A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
z = CSV.read("Real Data/JSE/Cleaned/JSECleanedTAQNPN.csv")
abs.(diff(log.(z[:,:MicroPrice])))
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
function AveAbsRet(data::DataFrame)
    # Extract the dates of the data
    dates = Date.(data[:,1])
    dates_unique = unique(dates)
    # Number of bins of bars each day
    nbins = length(findall(x->x==dates_unique[1], dates))-1
    # Create a matrix to store the normalised absolute intraday returns of each bin over the various days
    NormData = zeros(nbins, length(dates_unique))
    # Loop through each day of data
    for i in 1:length(dates_unique)
        # Extract the data for the day
        tempdata = data[findall(x->x==dates_unique[i], dates),:]
        # Get absolute intraday bar returns
        AbsDailyRet = abs.(diff(log.(tempdata[:,7])))
        # Get average absolute intraday bar returns
        AveAbsDailyRet = mean(AbsDailyRet)
        # Loop through each bin of the day
        for j in 1:nbins
            NormData[j, i] = AbsDailyRet[j] / AveAbsDailyRet
        end
    end
    # Initialise vector of aggregated absolute intraday returns
    AggregatedNormData = zeros(nbins, 1)
    # Compute the average absolute intraday returns for each bin and removing any empty bins
    for i in 1:nbins
        AggregatedNormData[i] = mean(filter(!isnan, NormData[i,:]))
    end
    times = Time.(data[2:(nbins+1),1])
    return times, AggregatedNormData
end

JSEBDT_AveAbsRet = AveAbsRet(JSETradesBars10minBDT)
JSEBDST_AveAbsRet = AveAbsRet(JSETradesBars10minBDST)

plot(JSEBDT_AveAbsRet[1], JSEBDT_AveAbsRet[2], seriestype = :bar, label = L"\textrm{BDT}", dpi = 300, alpha = 0.2, legend = :topright)
plot!(JSEBDST_AveAbsRet[1], JSEBDST_AveAbsRet[2], seriestype = :bar, label = L"\textrm{BDST}", alpha = 0.2)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Absolute Returns}")
# savefig("Assignment2/Plots/SBK_AveAbsRet.svg")

NPNBDT_AveAbsRet = AveAbsRet(NPNTradesBars10minBDT)
NPNBDST_AveAbsRet = AveAbsRet(NPNTradesBars10minBDST)

plot(NPNBDT_AveAbsRet[1], NPNBDT_AveAbsRet[2], seriestype = :bar, label = L"\textrm{BDT}", dpi = 300, alpha = 0.2, legend = :topright)
plot!(NPNBDST_AveAbsRet[1], NPNBDST_AveAbsRet[2], seriestype = :bar, label = L"\textrm{BDST}", alpha = 0.2)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Absolute Returns}")
# savefig("Assignment2/Plots/NPN_AveAbsRet.svg")

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
        if !isnan(line[10])
            # Check if current update is best bid or ask
            if line[2] == "BID" # If Bid
                # Get Spread = 2* (midprice - best bid)
                spread = 2 * abs(line[10] - line[3])
                # Push data into Master_df
                push!(Master_df, (line[1], spread))
            elseif line[2] == "ASK" # If ask
                # Get spread: 2 * (best ask - midprice)
                spread = 2 * abs(line[5] - line[10])
                # Push data into Master_df
                push!(Master_df, (line[1], spread, absReturn))
            end # We don't consider trades when making the spread
        end
    end
    return Master_df
end

## Function to compute the average spread normalised by
# the average daily spread
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
        # Only keep data within continuous trading
        start = DateTime(tempday) + Hour(9)
        close = DateTime(tempday) + Hour(16) + Minute(50)
        # Create the bars
        bars = collect(start:Minute(barsize):close)
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


SBKBDT_Spread = @pipe SBKBDT |> getSpread |> AveSpread(_, 10)
SBKBDST_Spread = @pipe SBKBDST |> getSpread |> AveSpread(_, 10)

plot(SBKBDT_Spread[1], SBKBDT_Spread[2], seriestype = :bar, label = L"\textrm{BDT}", dpi = 300, alpha = 0.2, legend = :topright)
plot!(SBKBDST_Spread[1], SBKBDST_Spread[2], seriestype = :bar, label = L"\textrm{BDST}", alpha = 0.2)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Spread}")
# savefig("Assignment2/Plots/SBK_Spread.svg")


NPNBDT_Spread = @pipe NPNBDT |> getSpread |> AveSpread(_, 10)
NPNBDST_Spread = @pipe NPNBDST |> getSpread |> AveSpread(_, 10)

plot(NPNBDT_Spread[1], NPNBDT_Spread[2], seriestype = :bar, label = L"\textrm{BDT}", dpi = 300, alpha = 0.2, legend = :topright)
plot!(NPNBDST_Spread[1], NPNBDST_Spread[2], seriestype = :bar, label = L"\textrm{BDST}", alpha = 0.2)
xlabel!(L"\textrm{Time of day}")
ylabel!(L"\textrm{Normalised Spread}")
# savefig("Assignment2/Plots/NPN_Spread.svg")
