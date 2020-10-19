## Author: Patrick Chang & Ivan Jerivich
# Script file to plot the price impact master curve of the 3 common equities
# from JSE and A2X combined together over the same period of 2019-01-02 to 2019-07-15.
# Looking for a possible universal master curve for the SA equity market (for common tickers)

#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim
using Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

# Tickers
tickers = ["SBK", "NPN", "SLM"]

# Read in the data
A2X_PriceImpact = load("Real Data/A2X/PriceImpact/A2X_PriceImpact.jld")
A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Real Data/JSE/PriceImpact/JSE_PriceImpact.jld")
JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]

#---------------------------------------------------------------------------
## Estimate γ and δ for SA market using 3 common equities from each exchange
# Note: δ = param[1]; γ = param[2]
# δ is for scaling of volume, γ is for scaling of impact
#---------------------------------------------------------------------------

function getErrorBuy(param; dataJSE = JSE_PriceImpact, dataA2X = A2X_PriceImpact, tickers = tickers, low = -1, up = 1)
    # Extract appropriate side
    ADV_JSE = dataJSE[3]; ADV_A2X = dataA2X[3]
    dataJSE = dataJSE[1]; dataA2X = dataA2X[1]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker in JSE
        for j in tickers
            # Pull out the data
            tempdata = dataJSE[j]
            C = ADV_JSE[j]
            # Find indicies for appropriate bin
            ind = findall(x-> x>xx[i-1] && x<=xx[i], tempdata[:,2])
            # Extract ω and Δp
            NormVol = tempdata[ind,2]
            Impact = tempdata[ind,1]
            # Remove any NaNs
            indofnan = findall(!isnan, Impact)
            NormVol = NormVol[indofnan]
            Impact = Impact[indofnan]
            # Rescale for by liquidity
            ω = [ω; mean(NormVol) / (C^δ)]
            Δp = [Δp; mean(Impact) .* (C^γ)]
        end
        # Loop through each ticker in A2X
        for j in tickers
            # Pull out the data
            tempdata = dataA2X[j]
            C = ADV_A2X[j]
            # Find indicies for appropriate bin
            ind = findall(x-> x>xx[i-1] && x<=xx[i], tempdata[:,2])
            # Extract ω and Δp
            NormVol = tempdata[ind,2]
            Impact = tempdata[ind,1]
            # Remove any NaNs
            indofnan = findall(!isnan, Impact)
            NormVol = NormVol[indofnan]
            Impact = Impact[indofnan]
            # Rescale for by liquidity
            ω = [ω; mean(NormVol) / (C^δ)]
            ω = filter(!isnan, ω)
            Δp = [Δp; mean(Impact) .* (C^γ)]
            Δp = filter(!isnan, Δp)
        end
        # Push components for the two-dimensional variance to be minimised
        push!(xbin, (std(ω) / mean(ω))^2)
        push!(ybin, (std(Δp) / mean(Δp))^2)
    end
    error = sum(xbin .+ ybin) / 20
    return error
end

function getErrorSell(param; dataJSE = JSE_PriceImpact, dataA2X = A2X_PriceImpact, tickers = tickers, low = -1, up = 1)
    # Extract appropriate side
    ADV_JSE = dataJSE[3]; ADV_A2X = dataA2X[3]
    dataJSE = dataJSE[2]; dataA2X = dataA2X[2]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker in JSE
        for j in tickers
            # Pull out the data
            tempdata = dataJSE[j]
            C = ADV_JSE[j]
            # Find indicies for appropriate bin
            ind = findall(x-> x>xx[i-1] && x<=xx[i], tempdata[:,2])
            # Extract ω and Δp
            NormVol = tempdata[ind,2]
            Impact = tempdata[ind,1]
            # Remove any NaNs
            indofnan = findall(!isnan, Impact)
            NormVol = NormVol[indofnan]
            Impact = Impact[indofnan]
            # Rescale for by liquidity
            ω = [ω; mean(NormVol) / (C^δ)]
            Δp = [Δp; mean(Impact) .* (C^γ)]
        end
        # Loop through each ticker in A2X
        for j in tickers
            # Pull out the data
            tempdata = dataA2X[j]
            C = ADV_A2X[j]
            # Find indicies for appropriate bin
            ind = findall(x-> x>xx[i-1] && x<=xx[i], tempdata[:,2])
            # Extract ω and Δp
            NormVol = tempdata[ind,2]
            Impact = tempdata[ind,1]
            # Remove any NaNs
            indofnan = findall(!isnan, Impact)
            NormVol = NormVol[indofnan]
            Impact = Impact[indofnan]
            # Rescale for by liquidity
            ω = [ω; mean(NormVol) / (C^δ)]
            ω = filter(!isnan, ω)
            Δp = [Δp; mean(Impact) .* (C^γ)]
            Δp = filter(!isnan, Δp)
        end
        # Push components for the two-dimensional variance to be minimised
        push!(xbin, (std(ω) / mean(ω))^2)
        push!(ybin, (std(Δp) / mean(Δp))^2)
    end
    error = sum(xbin .+ ybin) / 20
    return error
end

SACommonBuy = optimize(getErrorBuy, [0.3, 0.3])
SACommonBuyParam = SACommonBuy.minimizer

SACommonSell = optimize(getErrorSell, [0.3, 0.3])
SACommonSellParam = SACommonSell.minimizer

save("Computed Data/SACommonParams.jld", "SACommonBuyParam", SACommonBuyParam, "SACommonSellParam", SACommonSellParam)

## Load the parameters
#---------------------------------------------------------------------------

SACommonParams = load("Computed Data/SACommonParams.jld")
SACommonBuyParam = SACommonParams["SACommonBuyParam"]
SACommonSellParam = SACommonParams["SACommonSellParam"]

## Plot the results
#---------------------------------------------------------------------------

## Function to construct the plotting information for impact
function getMasterImpact(data::DataFrame, ADV::Float64, param::Vector; low = -1, up = 1)
    # Get parameters
    δ = param[1]; γ = param[2]
    data = data[findall(!isnan, data[:,1]),:]
    xx = 10 .^(range(low, up, length = 21))
    Δp = []
    ω = []
    for i in 2:length(xx)
        # Find indicies for appropriate bin
        ind = findall(x-> x>xx[i-1] && x<=xx[i], data[:,2])
        # Extract ω and Δp
        impacts = data[ind,1]
        normvol = data[ind,2]
        # Remove any NaNs
        indofnan = findall(!isnan, impacts)
        impacts = impacts[indofnan]
        normvol = normvol[indofnan]
        push!(Δp, mean(impacts) * (ADV^γ))
        push!(ω, mean(normvol) / (ADV^δ))
    end
    val_inds = setdiff(1:20, findall(iszero,Δp))
    val_inds = setdiff(val_inds, findall(isnan,Δp))
    return ω[val_inds], Δp[val_inds], val_inds
end

## Function to plot the price impact
function PlotMaster(param::Vector, side::Symbol; dataJSE = JSE_PriceImpact, dataA2X = A2X_PriceImpact, tickers = tickers, low = -1, up = 1)
    # Extract liquidity
    ADV_JSE = dataJSE[3]; ADV_A2X = dataA2X[3]
    # Extract appropriate side
    if side == :buy
        # Buyer-Initiated
        dataJSE = dataJSE[1]; dataA2X = dataA2X[1]
    elseif side == :sell
        # Seller-Initiated
        dataJSE = dataJSE[2]; dataA2X = dataA2X[2]
    end
    # Get Price Impact
    ImpactJSE = Dict()
    for i in tickers
        push!(ImpactJSE, i => getMasterImpact(dataJSE[i], ADV_JSE[i], param, low = low, up = up))
    end
    ImpactA2X = Dict()
    for i in tickers
        push!(ImpactA2X, i => getMasterImpact(dataA2X[i], ADV_A2X[i], param, low = low, up = up))
    end
    # Plot the values
    plot(ImpactJSE[tickers[1]][1], ImpactJSE[tickers[1]][2], marker = (4, 0.8), scale = :log10, dpi = 300,
    label = "JSE: SBK", legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400), color = :red)
    plot!(ImpactJSE[tickers[2]][1], ImpactJSE[tickers[2]][2], marker = (4, 0.8), scale = :log10, label = "JSE: NPN", color = :red, linestyle = :dash)
    plot!(ImpactJSE[tickers[3]][1], ImpactJSE[tickers[3]][2], marker = (4, 0.8), scale = :log10, label = "JSE: SLM", color = :red, linestyle = :dot)
    # for i in 2:length(tickers)
    #     plot!(ImpactJSE[tickers[i]][1], ImpactJSE[tickers[i]][2], marker = (4, 0.8), scale = :log10, label = "", color = :red)
    # end
    plot!(ImpactA2X[tickers[1]][1], ImpactA2X[tickers[1]][2], marker = (4, 0.8), scale = :log10, label = "A2X: SBK", color = :blue)
    plot!(ImpactA2X[tickers[2]][1], ImpactA2X[tickers[2]][2], marker = (4, 0.8), scale = :log10, label = "A2X: NPN", color = :blue, linestyle = :dash)
    plot!(ImpactA2X[tickers[3]][1], ImpactA2X[tickers[3]][2], marker = (4, 0.8), scale = :log10, label = "A2X: SLM", color = :blue, linestyle = :dot)
    # for i in 2:length(tickers)
    #     plot!(ImpactA2X[tickers[i]][1], ImpactA2X[tickers[i]][2], marker = (4, 0.8), scale = :log10, label = "", color = :blue)
    # end
    current()
    # Add appropriate label
    if side == :buy
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    end
end

PlotMaster(SACommonBuyParam, :buy)
# savefig("Plots/SACommonBuy.svg")

PlotMaster(SACommonSellParam, :sell)
# savefig("Plots/SACommonSell.svg")
