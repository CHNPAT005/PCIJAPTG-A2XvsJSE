### Title: Common impact
### Authors: Patrick Chang and Ivan Jericevich
### Function: Plot the price impact master curve of the common equities from JSE and A2X over the same period of 2019-01-02 to 2019-07-15
### Structure:
# 1. Preliminaries
# 2. Estimate γ and δ for using 3 common equities from each exchange
# 3. Visualization
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim, Statistics, LaTeXStrings, Distributions, StatsBase
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
ticker = ["SBK", "NPN", "SLM"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Estimate γ and δ for using 3 common equities from each exchange
function getJSEerrorBuy(param; data = JSE_PriceImpact, ticker = ticker, low = -1, up = 1)
    # Extract appropriate side
    ADV = data[3]
    data = data[1]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker
        for j in ticker
            # Pull out the data
            tempdata = data[j]
            C = ADV[j]
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
        # Push components for the two-dimensional variance to be minimised
        push!(xbin, (std(ω) / mean(ω))^2)
        push!(ybin, (std(Δp) / mean(Δp))^2)
    end
    error = sum(xbin .+ ybin) / 20
    return error
end
function getJSEerrorSell(param; data = JSE_PriceImpact, ticker = ticker, low = -1, up = 1)
    # Extract appropriate side
    ADV = data[3]
    data = data[2]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker
        for j in ticker
            # Pull out the data
            tempdata = data[j]
            C = ADV[j]
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
        # Push components for the two-dimensional variance to be minimised
        push!(xbin, (std(ω) / mean(ω))^2)
        push!(ybin, (std(Δp) / mean(Δp))^2)
    end
    error = sum(xbin .+ ybin) / 20
    return error
end
function getA2XerrorBuy(param; data = A2X_PriceImpact, ticker = ticker, low = -1, up = 1)
    # Extract appropriate side
    ADV = data[3]
    data = data[1]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker
        for j in ticker
            # Pull out the data
            tempdata = data[j]
            C = ADV[j]
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
function getA2XerrorSell(param; data = A2X_PriceImpact, ticker = ticker, low = -1, up = 1)
    # Extract appropriate side
    ADV = data[3]
    data = data[2]
    # Get parameters
    δ = param[1]; γ = param[2]
    # Initialise
    xbin = []; ybin = []
    xx = 10 .^(range(low, up, length = 21))
    # Loop through each bin
    for i in 2:length(xx)
        ω = []; Δp = []
        # Loop through each ticker
        for j in ticker
            # Pull out the data
            tempdata = data[j]
            C = ADV[j]
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
JSEBuy = optimize(getJSEerrorBuy, [0.3, 0.3]); JSEBuyParam = JSEBuy.minimizer
JSESell = optimize(getJSEerrorSell, [0.3, 0.3]); JSESellParam = JSESell.minimizer
A2XBuy = optimize(getA2XerrorBuy, [0.3, 0.3]); A2XBuyParam = A2XBuy.minimizer
A2XSell = optimize(getA2XerrorSell, [0.3, 0.3]); A2XSellParam = A2XSell.minimizer
save("Test Data/JSE/Price Impact/JSEParamsCommon.jld", "JSEBuyParam", JSEBuyParam, "JSESellParam", JSESellParam)
save("Test Data/A2X/Price Impact/A2XParamsCommon.jld", "A2XBuyParam", A2XBuyParam, "A2XSellParam", A2XSellParam)
JSEParams = load("Test Data/JSE/Price Impact/JSEParamsCommon.jld"); JSEBuyParam = JSEParams["JSEBuyParam"]; JSESellParam = JSEParams["JSESellParam"]
A2XParams = load("Test Data/A2X/Price Impact/A2XParamsCommon.jld"); A2XBuyParam = A2XParams["A2XBuyParam"]; A2XSellParam = A2XParams["A2XSellParam"]
#---------------------------------------------------------------------------


### 3. Visualization
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
    return ω[val_inds], Δp[val_inds]
end
function getPriceImpact(data::DataFrame; low = -3, up = 1)
    data = data[findall(!isnan, data[:,1]),:]
    xx = 10 .^(range(low, up, length = 21))
    Δp = []
    ω = []
    for i in 2:length(xx)
        ind = findall(x-> x>xx[i-1] && x<=xx[i], data[:,2])
        impacts = data[ind,1]
        normvol = data[ind,2]
        indofnan = findall(!isnan, impacts)
        impacts = impacts[indofnan]
        normvol = normvol[indofnan]
        push!(Δp, mean(impacts))
        push!(ω, mean(normvol))
    end
    val_inds = setdiff(1:length(Δp), findall(iszero,Δp))
    val_inds = setdiff(val_inds, findall(isnan,Δp))
    return ω[val_inds], Δp[val_inds]
end
function PlotMaster(data, ticker::Vector, param::Vector, side::Symbol; low = -1, up = 1)
    # Extract liquidity
    ADV = data[3]
    # Extract appropriate side
    if side == :buy
        # Buyer-Initiated
        data = data[1]
    elseif side == :sell
        # Seller-Initiated
        data = data[2]
    end
    # Get Price Impact
    ImpactRescaled = Dict()
    Impact = Dict()
    for i in ticker
        push!(ImpactRescaled, i => getMasterImpact(data[i], ADV[i], param, low = low, up = up))
        push!(Impact, i => getPriceImpact(data[i], low = low, up = up))
    end
    # Plot the values
    plot(ImpactRescaled[ticker[1]][1], ImpactRescaled[ticker[1]][2], marker = (4, 0.8), scale = :log10, color = :blue, dpi = 300, label = ticker[1], legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400))
    plot!(Impact[ticker[1]][1], Impact[ticker[1]][2], marker = (4, :white, stroke(1, :blue)), scale = :log10, color = :blue, label = "", line = (:dot, 1))
    plot!(ImpactRescaled[ticker[2]][1], ImpactRescaled[ticker[2]][2], marker = (4, 0.8), scale = :log10, label = ticker[2], color= :red)
    plot!(Impact[ticker[2]][1], Impact[ticker[2]][2], marker = (4, :white, stroke(1, :red)), scale = :log10, color = :red, label = "", line = (:dot, 1))
    plot!(ImpactRescaled[ticker[3]][1], ImpactRescaled[ticker[3]][2], marker = (4, 0.8), scale = :log10, label = ticker[3], color= :green)
    plot!(Impact[ticker[3]][1], Impact[ticker[3]][2], marker = (4, :white, stroke(1, :green)), scale = :log10, color = :green, label = "", line = (:dot, 1))
    # Add appropriate label
    if side == :buy
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    end
end
PlotMaster(JSE_PriceImpact, ticker, JSEBuyParam, :buy); savefig("Figures/JSECommonBuy.pdf")
PlotMaster(JSE_PriceImpact, ticker, JSESellParam, :sell); savefig("Figures/JSECommonSell.pdf")
PlotMaster(A2X_PriceImpact, ticker, A2XBuyParam, :buy); savefig("Figures/A2XCommonBuy.pdf")
PlotMaster(A2X_PriceImpact, ticker, A2XSellParam, :sell); savefig("Figures/A2XCommonSell.pdf")
#---------------------------------------------------------------------------
