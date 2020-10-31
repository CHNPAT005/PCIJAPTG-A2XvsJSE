### Title: Order-book properties
### Authors: Patrick Chang and Ivan Jericevich
### Function: Investigate order flow auto-correlation and micro-price return auto-correlation (removing overnight returns)
### Structure:
# 1. Preliminaries
# 2. Order-flow auto-correlation
# 3. Micro-price return auto-correlation
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings, Pipe, StatsBase, Distributions
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
A2X = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_NPN.csv"); JSE = CSV.read("Test Data/JSE/Clean/JSECleanedTAQNPN.csv")
#---------------------------------------------------------------------------


### 2. Order-flow auto-correlation
function getA2XTradeSigns(data::DataFrame) # Obtain the Lee-Ready classification for A2X
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
function getJSETradeSigns(data::DataFrame) # Obtain the Lee-Ready classification for JSE
    tradeinds = findall(x -> x == "TRADE", data[:,2])
    inferredclassification = data[tradeinds,12]
    return inferredclassification
end
function ConvertTradeSigns(data::Vector) # Take in vector of strings with trade classification as buyer- or seller-initiated and convert the strings to +1 and -1.
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
A2XOrderFlow = @pipe A2X |> getA2XTradeSigns |> ConvertTradeSigns
JSEOrderFlow = @pipe JSE |> getJSETradeSigns |> convert(Vector, _) |> ConvertTradeSigns
lags = 1000
A2X_OF = plot(1:lags, autocor(A2XOrderFlow, 1:lags), seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black)
hline!(A2X_OF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(A2XOrderFlow))], color = :blue)
hline!(A2X_OF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(A2XOrderFlow))], color = :blue)
plot!(A2X_OF, 1:lags, autocor(A2XOrderFlow, 1:lags), xscale = :log10,; inset = (1, bbox(0.58,0.0,0.4,0.4)), subplot = 2, label = "", xlabel = L"\textrm{Lag} (\log_{10})", ylabel = L"\textrm{ACF}", guidefontsize = 8, color = :black); savefig(A2X_OF, "Figures/A2X_OF.pdf")
JSE_OF = plot(1:lags, autocor(JSEOrderFlow, 1:lags), seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black)
hline!(JSE_OF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(JSEOrderFlow))], color = :red)
hline!(JSE_OF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(JSEOrderFlow))], color = :red)
plot!(JSE_OF, 1:lags, autocor(JSEOrderFlow, 1:lags), xscale = :log10, inset = (1, bbox(0.58,0.0,0.4,0.4)), subplot = 2, label = "", xlabel = L"\textrm{Lag} (\log_{10})", ylabel = L"\textrm{ACF}", guidefontsize = 8, color = :black); savefig(JSE_OF, "Figures/JSE_OF.pdf")
#---------------------------------------------------------------------------


### 3. Micro-price return auto-correlation
function getA2XReturns(data::DataFrame) # Extract the micro-price returns from A2X
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
function getJSEReturns(data::DataFrame) # Extract the micro-price returns from JSE
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
A2XRets = getA2XReturns(A2X); JSERets = getJSEReturns(JSE) # Compute the micro-price returns
lags = 100
A2X_Rets_ACF = plot(1:lags, autocor(convert.(Float64, A2XRets), 1:lags), seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black, label = "")
hline!(A2X_Rets_ACF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(A2XRets))], color = :blue, label = "")
hline!(A2X_Rets_ACF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(A2XRets))], color = :blue, label = ""); savefig(A2X_Rets_ACF, "Figures/A2X_Rets_ACF.pdf")
JSE_Rets_ACF = plot(1:lags, autocor(convert.(Float64, JSERets), 1:lags), seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", dpi = 300, color = :black, label = "")
hline!(JSE_Rets_ACF, [quantile(Normal(), (1+0.95)/2) / sqrt(length(JSERets))], color = :red, label = "")
hline!(JSE_Rets_ACF, [quantile(Normal(), (1-0.95)/2) / sqrt(length(JSERets))], color = :red, label = ""); savefig(JSE_Rets_ACF, "Figures/JSE_Rets_ACF.pdf")
#---------------------------------------------------------------------------
