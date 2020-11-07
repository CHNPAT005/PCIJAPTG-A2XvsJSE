### Title: Price impact
### Authors: Patrick Chang and Ivan Jericevich
### Function: Visualise the price impact curves of JSE and A2X equities over the same period
### Structure:
# 1. Preliminaries
# 2. Compute price impact curves
# 3. Visualization
#---------------------------------------------------------------------------


### 1. Preamble
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings, StatsBase, ColorSchemes
cd("C:/Users/Ivan/Documents/PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Compute price impact curves
function getPriceImpact(data::DataFrame; low = -3, up = 1)
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
        push!(Δp, mean(impacts))
        push!(ω, mean(normvol))
    end
    val_inds = setdiff(1:20, findall(iszero,Δp))
    val_inds = setdiff(val_inds, findall(isnan,Δp))
    return ω[val_inds], Δp[val_inds], val_inds
end
function getBootImpact(data::DataFrame; low = -3, up = 1)
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
        push!(Δp, mean(impacts))
        push!(ω, mean(normvol))
    end
    val_inds = setdiff(1:20, findall(iszero,Δp))
    val_inds = setdiff(val_inds, findall(isnan,Δp))
    return ω[val_inds], Δp[val_inds], val_inds
end
function PlotBootstrap(data, M::Int, ticker::Vector, side::Symbol, cutoff::Float64; low = -3, up = 1)
    # Extract appropriate side
    if side == :buy
        # Buyer-Initiated
        data = data[1]
    elseif side == :sell
        # Seller-Initiated
        data = data[2]
    end
    # Get Price Impact
    realImpact = Dict()
    for i in ticker
        push!(realImpact, i => getPriceImpact(data[i], low = low, up = up))
    end
    # Initialize re-scaled impact for the bootstrap iteration
    tempmasterΔp = fill(NaN, 10, 20, M)
    tempmasterω = fill(NaN, 10, 20, M)
    # Loop through the M bootstraps
    @showprogress "Computing:" for i in 1:M
        # Get Price Impact
        Impact = Dict()
        for j in ticker
            push!(Impact, j => getBootImpact(data[j], low = low, up = up))
        end
        # Populate the re-scaled impact for each bootstrap iteration in a 3D matrix
        for j in 1:length(ticker)
            tempmasterΔp[j, Impact[ticker[j]][3], i] = Impact[ticker[j]][2]
            tempmasterω[j, Impact[ticker[j]][3], i] = Impact[ticker[j]][1]
        end
    end
    # return tempmasterΔp
    Boots = Dict()
    for j in 1:size(tempmasterΔp, 1)
        temp = []
        tempMean = []; tempStd = []
        for k in 1:size(tempmasterΔp, 2)
            # Get the non-nan indeces on the third dimension (all the bootstrap samples for a security and volume bin)
            nonnanindeces = findall(!isnan, tempmasterΔp[j, k, :])
            ave = mean(tempmasterΔp[j, k, nonnanindeces])
            dev = std(tempmasterΔp[j, k, nonnanindeces])
            if !isnan(ave) && !isnan(dev)
                push!(tempMean, ave)
                push!(tempStd, dev)
            end

        end
        push!(temp, tempMean); push!(temp, tempStd)
        push!(Boots, ticker[j] => temp)
    end
    # Plot the values
    plot(realImpact[ticker[1]][1], realImpact[ticker[1]][2], scale = :log10, dpi = 300, label = "", legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400), fillrange = max.(realImpact[ticker[1]][2] .- 1.96 .* Boots[ticker[1]][2], cutoff), fillalpha = 0.2, fillcolor = 1, palette = ColorSchemes.tab10.colors)
    plot!(realImpact[ticker[1]][1], realImpact[ticker[1]][2], scale = :log10, label = "", legend = :outertopright, fillrange = max.(realImpact[ticker[1]][2] .+ 1.96 .* Boots[ticker[1]][2], cutoff), fillalpha = 0.2, fillcolor = 1)
    plot!(realImpact[ticker[1]][1], realImpact[ticker[1]][2], marker = (4, 0.8), scale = :log10, label = ticker[1], legend = :outertopright, markercolor = 1, markerstrokecolor = 1, linecolor = 1)
    for i in 2:length(ticker)
        plot!(realImpact[ticker[i]][1], realImpact[ticker[i]][2], scale = :log10, dpi = 300, label = "", legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700,400), fillrange = max.(realImpact[ticker[i]][2] .- 1.96 .* Boots[ticker[i]][2], cutoff), fillalpha = 0.2, fillcolor = i)
        plot!(realImpact[ticker[i]][1], realImpact[ticker[i]][2], scale = :log10, label = "", legend = :outertopright, fillrange = max.(realImpact[ticker[i]][2] .+ 1.96 .* Boots[ticker[i]][2], cutoff), fillalpha = 0.2, fillcolor = i)
        plot!(realImpact[ticker[i]][1], realImpact[ticker[i]][2], marker = (4, 0.8), scale = :log10, label = ticker[i], legend = :outertopright, markercolor = i, markerstrokecolor = i, linecolor = i)
    end
    current()
    # Add appropriate label
    if side == :buy
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^*")
        ylabel!(L"\Delta p^*")
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^*")
        ylabel!(L"\Delta p^*")
    end
end
#---------------------------------------------------------------------------


### 3. Visualization
PlotBootstrap(JSE_PriceImpact, 1000, JSE_tickers, :buy, 10^(-6))
savefig("Figures/JSEImpactBuy.svg")

PlotBootstrap(JSE_PriceImpact, 1000, JSE_tickers, :sell, 10^(-6))
savefig("Figures/JSEImpactSell.svg")

PlotBootstrap(A2X_PriceImpact, 1000, A2X_tickers, :buy, 10^(-5))
savefig("Figures/A2XImpactBuy.svg")

PlotBootstrap(A2X_PriceImpact, 1000, A2X_tickers, :sell, 10^(-5))
savefig("Figures/A2XImpactSell.svg")
#---------------------------------------------------------------------------
