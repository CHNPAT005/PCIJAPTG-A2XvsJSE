## Author: Patrick Chang & Ivan Jerivich
# Script file to plot the price impact master curve of 10 equities
# from JSE and A2X over the same period of 2019-01-02 to 2019-07-15.

#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim
using Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

# JSE tickers
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]

# A2X tickers
A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]

# Read in the data
A2X_PriceImpact = load("Real Data/A2X/PriceImpact/A2X_PriceImpact.jld")
A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Real Data/JSE/PriceImpact/JSE_PriceImpact.jld")
JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]

#---------------------------------------------------------------------------
## Estimate γ and δ for each exchange
# Note: δ = param[1]; γ = param[2]
# δ is for scaling of volume, γ is for scaling of impact
#---------------------------------------------------------------------------

function getJSEerrorBuy(param; data = JSE_PriceImpact, ticker = JSE_tickers, low = -1, up = 1)
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

function getJSEerrorSell(param; data = JSE_PriceImpact, ticker = JSE_tickers, low = -1, up = 1)
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

function getA2XerrorBuy(param; data = A2X_PriceImpact, ticker = A2X_tickers, low = -1, up = 1)
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

function getA2XerrorSell(param; data = A2X_PriceImpact, ticker = A2X_tickers, low = -1, up = 1)
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


JSEBuy = optimize(getJSEerrorBuy, [0.3, 0.3])
JSEBuyParam = JSEBuy.minimizer

JSESell = optimize(getJSEerrorSell, [0.3, 0.3])
JSESellParam = JSESell.minimizer

A2XBuy = optimize(getA2XerrorBuy, [0.3, 0.3])
A2XBuyParam = A2XBuy.minimizer

A2XSell = optimize(getA2XerrorSell, [0.3, 0.3])
A2XSellParam = A2XSell.minimizer

save("Computed Data/JSEParams.jld", "JSEBuyParam", JSEBuyParam, "JSESellParam", JSESellParam)
save("Computed Data/A2XParams.jld", "A2XBuyParam", A2XBuyParam, "A2XSellParam", A2XSellParam)

## Load the parameters
#---------------------------------------------------------------------------

JSEParams = load("Computed Data/JSEParams.jld")
JSEBuyParam = JSEParams["JSEBuyParam"]
JSESellParam = JSEParams["JSESellParam"]
A2XParams = load("Computed Data/A2XParams.jld")
A2XBuyParam = A2XParams["A2XBuyParam"]
A2XSellParam = A2XParams["A2XSellParam"]

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
    return ω[val_inds], Δp[val_inds]
end

## Function to plot the price impact
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
    Impact = Dict()
    for i in ticker
        push!(Impact, i => getMasterImpact(data[i], ADV[i], param, low = low, up = up))
    end
    # Plot the values
    plot(Impact[ticker[1]][1], Impact[ticker[1]][2], marker = (4, 0.8), scale = :log10, dpi = 300,
    label = ticker[1], legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400))
    for i in 2:length(ticker)
        plot!(Impact[ticker[i]][1], Impact[ticker[i]][2], marker = (4, 0.8), scale = :log10, label = ticker[i])
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

## Obtain results

PlotMaster(JSE_PriceImpact, JSE_tickers, JSEBuyParam, :buy)
# savefig("Plots/JSEMasterBuy.svg")

PlotMaster(JSE_PriceImpact, JSE_tickers, JSESellParam, :sell)
# savefig("Plots/JSEMasterSell.svg")

PlotMaster(A2X_PriceImpact, A2X_tickers, A2XBuyParam, :buy)
# savefig("Plots/A2XMasterBuy.svg")

PlotMaster(A2X_PriceImpact, A2X_tickers, A2XSellParam, :sell)
# savefig("Plots/A2XMasterSell.svg")
