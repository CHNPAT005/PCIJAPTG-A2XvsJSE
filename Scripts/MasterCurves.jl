### Title: Master curves
### Authors: Patrick Chang and Ivan Jericevich
### Function: generate master curves for 10 equities in each exchange
### Structure:
# 1. Preliminaries
# 2. Estimate γ and δ for each exchange
# 3. Visualization
# 4. Compare the master curves from the exchanges
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim, Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase, Pipe
cd("C:/Users/Ivan/Documents/PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
# A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
# JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
A2X_PriceImpact = load("Real Data/A2X/PriceImpact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Real Data/JSE/PriceImpact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Estimate γ and δ for each exchange and for buyer and seller-initiated trades
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
JSEBuy = optimize(getJSEerrorBuy, [0.3, 0.3]); JSEBuyParam = JSEBuy.minimizer
JSESell = optimize(getJSEerrorSell, [0.3, 0.3]); JSESellParam = JSESell.minimizer
A2XBuy = optimize(getA2XerrorBuy, [0.3, 0.3]); A2XBuyParam = A2XBuy.minimizer
A2XSell = optimize(getA2XerrorSell, [0.3, 0.3]); A2XSellParam = A2XSell.minimizer
save("Computed Data/JSEParams.jld", "JSEBuyParam", JSEBuyParam, "JSESellParam", JSESellParam)
save("Computed Data/A2XParams.jld", "A2XBuyParam", A2XBuyParam, "A2XSellParam", A2XSellParam)
JSEParams = load("Computed Data/JSEParams.jld"); JSEBuyParam = JSEParams["JSEBuyParam"]; JSESellParam = JSEParams["JSESellParam"]
A2XParams = load("Computed Data/A2XParams.jld"); A2XBuyParam = A2XParams["A2XBuyParam"]; A2XSellParam = A2XParams["A2XSellParam"]
#---------------------------------------------------------------------------

### 3. Visualization
# Using only raw data
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
    # plot(Impact[ticker[1]][1], Impact[ticker[1]][2], marker = (4, 0.8), scale = :log10, dpi = 300, label = ticker[1], legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400))
    # plot!()
    for i in 1:length(ticker)
        plot!(Impact[ticker[i]][1], Impact[ticker[i]][2], marker = (4, 0.8), scale = :log10, label = ticker[i], palette = ColorSchemes.Dark2_5.colors)
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

# Bootstrapping raw data
function getBootImpact(data::DataFrame, ADV::Float64, param::Vector; low = -1, up = 1)
    # Get parameters
    δ = param[1]; γ = param[2]
    data = data[findall(!isnan, data[:,1]),:]
    N = size(data, 1)
    bootsamples = sample(1:N, N, replace = true)
    data = data[bootsamples,:]
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

function PlotBootstrap(data, M::Int, ticker::Vector, param::Vector, side::Symbol, col::Symbol; low = -1, up = 1)
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
    # Initialize bootstrap master curves
    masterscurvesΔp = fill(NaN, M, 20)
    masterscurvesω = fill(NaN, M, 20)
    # Loop through the M bootstraps
    @showprogress "Computing:" for i in 1:M
        # Get Price Impact
        Impact = Dict()
        for j in ticker
            push!(Impact, j => getBootImpact(data[j], ADV[j], param, low = low, up = up))
        end
        # Initialize re-scaled impact for the bootstrap iteration
        tempmasterΔp = fill(NaN, 10, 20)
        tempmasterω = fill(NaN, 10, 20)
        # Populate the re-scaled impact for the bootstrap iteration
        for j in 1:length(ticker)
            tempmasterΔp[j, Impact[ticker[j]][3]] = Impact[ticker[j]][2]
            tempmasterω[j, Impact[ticker[j]][3]] = Impact[ticker[j]][1]
        end
        # Populate the bootstrap master curves
        for j in 1:20
            masterscurvesΔp[i,j] = @pipe tempmasterΔp[:,j] |> filter(!isnan,_) |> mean
            masterscurvesω[i,j] = @pipe tempmasterω[:,j] |> filter(!isnan,_) |> mean
        end
    end
    # Plot the values
    plot(mean(masterscurvesω, dims=1)', mean(masterscurvesΔp, dims=1)', ribbon = 1.96 .* std(masterscurvesΔp, dims=1)', color = col, fillalpha = 0.3, scale = :log10, dpi = 300, label = "Average", legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400))
    # Add appropriate label
    if side == :buy
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^* / C^{\delta}")
        ylabel!(L"\Delta p^* C^{\gamma}")
    end
end

PlotBootstrap(A2X_PriceImpact, 1000, A2X_tickers, A2XBuyParam, :buy, :blue)
PlotMaster(A2X_PriceImpact, A2X_tickers, A2XBuyParam, :buy)
savefig("Figures/A2XMasterBuy.svg")

PlotBootstrap(A2X_PriceImpact, 1000, A2X_tickers, A2XSellParam, :sell, :blue)
PlotMaster(A2X_PriceImpact, A2X_tickers, A2XSellParam, :sell)
savefig("Figures/A2XMasterSell.svg")

PlotBootstrap(JSE_PriceImpact, 1000, JSE_tickers, JSEBuyParam, :buy, :red)
PlotMaster(JSE_PriceImpact, JSE_tickers, JSEBuyParam, :buy)
savefig("Figures/JSEMasterBuy.svg")

PlotBootstrap(JSE_PriceImpact, 1000, JSE_tickers, JSESellParam, :sell, :red)
PlotMaster(JSE_PriceImpact, JSE_tickers, JSESellParam, :sell)
savefig("Figures/JSEMasterSell.svg")
