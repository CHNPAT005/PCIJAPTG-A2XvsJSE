## Author: Patrick Chang and Ivan Jerivich
# Script file to investigate:
# 1) Order flow auto-correlation on A2X and JSE
#       Note that the order-flow uses the Lee-Ready classification
# 2) Micro-price return auto-correlation on A2X and JSE
#       Note that we compute returns as price fluctuations
#       Also we remove overnight returns

## Preamble

using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings, Pipe, StatsBase, Distributions

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

A2X = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")
JSE = CSV.read("Real Data/JSE/Cleaned/JSECleanedTAQNPN.csv")

# Component 1. Order-flow auto-correlation.
#---------------------------------------------------------------------------
# Function to obtain the Lee-Ready classification for A2X
function getA2XTradeSigns(data::DataFrame)
    # Initialise inferred classification
    inferredclassification = String[]
    # Split the data into individual days
    dates = Date.(data[:,2])
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01") && x <= Date("2019-07-15"), dates_unique)
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,3])
        # Extract trade values
        TradeValues = tempdata[tradeinds, 8]
        for j in 1:size(TradeValues)[1]
            # Get trade value
            TradeValue = tempdata[tradeinds[j], 8]
            MidQuoteBeforeTrade = tempdata[findlast(!isnan, tempdata[1:tradeinds[j], 12]), 12]
            # Perform the logic checks: begin with quote rule
            if TradeValue > MidQuoteBeforeTrade
                # Transaction is higher than mid price => BuyerInitiated
                push!(inferredclassification, "BuyerInitiated")
            elseif TradeValue < MidQuoteBeforeTrade
                # Transaction is lower than mid price => SellerInitiated
                push!(inferredclassification, "SellerInitiated")
            elseif TradeValue == MidQuoteBeforeTrade
                # Quote rule failed, go to tick rule
                if j > 1
                    if TradeValues[j] > TradeValues[j-1]
                        # If trade is higher than previous trade => BuyerInitiated
                        push!(inferredclassification, "BuyerInitiated")
                    elseif TradeValues[j] < TradeValues[j-1]
                        # If trade is lower than previous trade => SellerInitiated
                        push!(inferredclassification, "SellerInitiated")

                    elseif TradeValues[j] == TradeValues[j-1]
                        # If trades are the same, find last trade that was different
                        indTradeLast = findlast(x -> x != TradeValues[j], TradeValues[1:j])
                        if isnothing(indTradeLast)
                            # No classification, all trades before was the same
                            inferredclassification[j] = ""
                        else
                            # Compare against last trade that was different
                            TradeLast = TradeValues[indTradeLast]
                            if TradeValues[j] > TradeLast
                                # If trade is higher than previous trade => BuyerInitiated
                                inferredclassification[j] = "BuyerInitiated"
                            elseif TradeValues[j] < TradeLast
                                # If trade is lower than previous trade => SellerInitiated
                                inferredclassification[j] = "SellerInitiated"
                            end
                        end
                    end
                else
                    # If first trade can't be classified
                    push!(inferredclassification, "")
                end
            end
        end
    end
    return inferredclassification
end

# Function to obtain the Lee-Ready classification for JSE
# (already included in cleaned data just need to extract it)
function getJSETradeSigns(data::DataFrame)
    tradeinds = findall(x -> x == "TRADE", data[:,2])
    inferredclassification = data[tradeinds,12]
    return inferredclassification
end

# Function to take in vector of strings with trade classification as buyer- or seller-initiated
# and converting the strings to +1 and -1.
# Note that BuyerInitiated = +1, SellerInitiated = -1
function ConvertTradeSigns(data::Vector)
    # Initialise vector
    orderflow = Float64[]
    for i in 1:length(data)
        # Loop through the data
        if data[i] == "BuyerInitiated"
            # BuyerInitiated => 1
            push!(orderflow, 1)
        elseif data[i] == "SellerInitiated"
            # SellerInitiated => -1
            push!(orderflow, -1)
        end
        # if not classification then do not append a sign
    end
    return orderflow
end

## Compute the Order Flow

A2XOrderFlow = @pipe A2X |> getA2XTradeSigns |> ConvertTradeSigns

JSEOrderFlow = @pipe JSE |> getJSETradeSigns |> convert(Vector, _) |> ConvertTradeSigns

## Visualise the Order Flow autocorrelation
lags = 1000

A2X_OF = plot(1:lags, autocor(A2XOrderFlow, 1:lags),
seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black)
hline!(A2X_OF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(A2XOrderFlow))], color = :blue)
hline!(A2X_OF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(A2XOrderFlow))], color = :blue)
plot!(A2X_OF, 1:lags, autocor(A2XOrderFlow, 1:lags), xscale = :log10,
inset = (1, bbox(0.58,0.0,0.4,0.4)), subplot = 2, label = "",
xlabel = L"\textrm{Lag} (\log_{10})", ylabel = L"\textrm{ACF}", guidefontsize = 8, color = :black)
# savefig(A2X_OF, "Plots/A2X_OF.svg")

JSE_OF = plot(1:lags, autocor(JSEOrderFlow, 1:lags),
seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black)
hline!(JSE_OF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(JSEOrderFlow))], color = :red)
hline!(JSE_OF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(JSEOrderFlow))], color = :red)
plot!(JSE_OF, 1:lags, autocor(JSEOrderFlow, 1:lags), xscale = :log10,
inset = (1, bbox(0.58,0.0,0.4,0.4)), subplot = 2, label = "",
xlabel = L"\textrm{Lag} (\log_{10})", ylabel = L"\textrm{ACF}", guidefontsize = 8, color = :black)
# savefig(JSE_OF, "Plots/JSE_OF.svg")

# Component 2. Micro-price return autocorrelation
#---------------------------------------------------------------------------
# Note that here the returns are computed as the price fluctuations
# i.e. r(t_{k}) = log(p(t_{k+1})) - log(p(t_{k})).
# Also note that overnight returns are removed

# Function to extract the micro-price returns from A2X
function getA2XReturns(data::DataFrame)
    # Initialise returns
    returns = String[]
    # Split the data into individual days
    dates = Date.(data[:,2])
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01") && x <= Date("2019-07-15"), dates_unique)
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Extract the micro-prices
        MP = tempdata[:,11]
        # Compute and append returns
        returns = [returns; diff(log.(filter(!isnan, MP)))]
    end
    return returns
end

# Function to extract the micro-price returns from JSE
function getJSEReturns(data::DataFrame)
    # Initialise returns
    returns = String[]
    # Split the data into individual days
    dates = Date.(data[:,1])
    dates_unique = unique(dates)
    dates_unique = filter(x -> x >= Date("2019-01-01") && x <= Date("2019-07-15"), dates_unique)
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Extract the micro-prices
        MP = tempdata[:,9]
        # Compute and append returns
        returns = [returns; diff(log.(filter(!isnan, MP)))]
    end
    return returns
end

## Compute the micro-price returns

A2XRets = getA2XReturns(A2X)

JSERets = getJSEReturns(JSE)

## Visualise the micro-price return autocorrelation

lags = 100
A2X_Rets_ACF = plot(1:lags, autocor(convert.(Float64, A2XRets), 1:lags),
seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black, label = "")
hline!(A2X_Rets_ACF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(A2XRets))], color = :blue, label = "")
hline!(A2X_Rets_ACF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(A2XRets))], color = :blue, label = "")
# savefig(A2X_Rets_ACF, "Plots/A2X_Rets_ACF.svg")

JSE_Rets_ACF = plot(1:lags, autocor(convert.(Float64, JSERets), 1:lags),
seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black, label = "")
hline!(JSE_Rets_ACF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(JSERets))], color = :red, label = "")
hline!(JSE_Rets_ACF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(JSERets))], color = :red, label = "")
# savefig(JSE_Rets_ACF, "Plots/JSE_Rets_ACF.svg")
