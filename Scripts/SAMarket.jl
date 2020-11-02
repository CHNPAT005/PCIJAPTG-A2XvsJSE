### Title: Master curves for 10 equities
### Authors: Patrick Chang and Ivan Jericevich
### Function: Looking for a possible universal master curve for the SA equity market
### Structure:
# 1. Preliminaries
# 2. Estimate γ and δ for SA market using 10 equities from each exchange
# 3. Visualisation
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim, Statistics, LaTeXStrings, StatsBase
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Estimate γ and δ for SA market using 10 equities from each exchange
function getErrorBuy(param; dataJSE = JSE_PriceImpact, JSEticker = JSE_tickers, dataA2X = A2X_PriceImpact, A2Xticker = A2X_tickers, low = -1, up = 1) # δ is for scaling of volume, γ is for scaling of impact
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
        for j in JSEticker
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
        for j in A2Xticker
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
function getErrorSell(param; dataJSE = JSE_PriceImpact, JSEticker = JSE_tickers, dataA2X = A2X_PriceImpact, A2Xticker = A2X_tickers, low = -1, up = 1)
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
        for j in JSEticker
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
        for j in A2Xticker
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
SABuy = optimize(getErrorBuy, [0.3, 0.3]); SABuyParam = SABuy.minimizer
SASell = optimize(getErrorSell, [0.3, 0.3]); SASellParam = SASell.minimizer
save("Computed Data/SAParams.jld", "SABuyParam", SABuyParam, "SASellParam", SASellParam)
SAParams = load("Computed Data/SAParams.jld"); SABuyParam = SAParams["SABuyParam"]; SASellParam = SAParams["SASellParam"]
#---------------------------------------------------------------------------


### 3. Visualisation
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
function PlotMaster(param::Vector, side::Symbol; dataJSE = JSE_PriceImpact, JSEticker = JSE_tickers, dataA2X = A2X_PriceImpact, A2Xticker = A2X_tickers, low = -1, up = 1)
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
    for i in JSEticker
        push!(ImpactJSE, i => getMasterImpact(dataJSE[i], ADV_JSE[i], param, low = low, up = up))
    end
    ImpactA2X = Dict()
    for i in A2Xticker
        push!(ImpactA2X, i => getMasterImpact(dataA2X[i], ADV_A2X[i], param, low = low, up = up))
    end
    # Plot the values
    plot(ImpactJSE[JSEticker[1]][1], ImpactJSE[JSEticker[1]][2], marker = (4, 0.8), scale = :log10, dpi = 300,
    label = "JSE", legend = :outertopright, legendtitle = L"\textrm{Exchange}", size = (700,400), color = :red)
    for i in 2:length(JSEticker)
        plot!(ImpactJSE[JSEticker[i]][1], ImpactJSE[JSEticker[i]][2], marker = (4, 0.8), scale = :log10, label = "", color = :red)
    end
    plot!(ImpactA2X[A2Xticker[1]][1], ImpactA2X[A2Xticker[1]][2], marker = (4, 0.8), scale = :log10, label = "A2X", color = :blue)
    for i in 2:length(A2Xticker)
        plot!(ImpactA2X[A2Xticker[i]][1], ImpactA2X[A2Xticker[i]][2], marker = (4, 0.8), scale = :log10, label = "", color = :blue)
    end
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
PlotMaster(SABuyParam, :buy); savefig("Figures/SABuy.svg")
PlotMaster(SASellParam, :sell); savefig("Figures/SASell.svg")
#---------------------------------------------------------------------------
