### Title: Data Cleaning - A2X
### Authors: Patrick Chang and Ivan Jericevich
### Function: Plot the price impact (with direct cost incorporated) for 10 equities from JSE and A2X over the same period of 2019-01-02 to 2019-07-15
### Structure:
# 1. Preliminaries
# 2. Compute and plot costs
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Compute and plot costs
function getCosts(data::DataFrame, type::Symbol; low = -3, up = 1)
    data = data[findall(!isnan, data[:,1]),:]
    xx = 10 .^(range(low, up, length = 21))
    Δp = []
    ω = []
    for i in 2:length(xx)
        ind = findall(x-> x>xx[i-1] && x<=xx[i], data[:,2])
        costs = data[ind,type]
        normvol = data[ind,2]
        indofnan = findall(!isnan, costs)
        costs = costs[indofnan]
        normvol = normvol[indofnan]
        push!(Δp, mean(costs))
        push!(ω, mean(normvol))
    end
    val_inds = setdiff(1:length(Δp), findall(iszero,Δp))
    val_inds = setdiff(val_inds, findall(isnan,Δp))
    return ω[val_inds], Δp[val_inds]
end
function PlotCosts(type::Symbol, side::Symbol; dataJSE = JSE_PriceImpact, JSEticker = JSE_tickers, dataA2X = A2X_PriceImpact, A2Xticker = A2X_tickers, low = -1, up = 1)
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
        push!(ImpactJSE, i => getCosts(dataJSE[i], type, low = low, up = up))
    end
    ImpactA2X = Dict()
    for i in A2Xticker
        push!(ImpactA2X, i => getCosts(dataA2X[i], type, low = low, up = up))
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
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^*")
        if type == :Cost
            ylabel!(L"\Delta c^*")
        elseif type == :DirectCost
            ylabel!(L"\Delta \textrm{DC}^*")
        elseif type == :Spread
            ylabel!(L"\Delta s^*")
        elseif type == :Impact
            ylabel!(L"\Delta p^*")
        end
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^*")
        if type == :Cost
            ylabel!(L"\Delta c^*")
        elseif type == :DirectCost
            ylabel!(L"\Delta \textrm{DC}^*")
        elseif type == :Spread
            ylabel!(L"\Delta s^*")
        elseif type == :Impact
            ylabel!(L"\Delta p^*")
        end
    end
end
PlotCosts(:Cost, :buy); savefig("Figures/BuyCost.svg")
PlotCosts(:Cost, :sell); savefig("Figures/SellCost.svg")
PlotCosts(:DirectCost, :buy); savefig("Figures/DirectCostBuy.svg")
PlotCosts(:DirectCost, :sell); savefig("Figures/DirectCostSell.svg")
PlotCosts(:Spread, :buy); savefig("Figures/SpreadCostBuy.svg")
PlotCosts(:Spread, :sell); savefig("Figures/SpreadCostSell.svg")
PlotCosts(:Impact, :buy); savefig("Figures/ImpactCostBuy.svg")
PlotCosts(:Impact, :sell); savefig("Figures/ImpactCostSell.svg")
#---------------------------------------------------------------------------
