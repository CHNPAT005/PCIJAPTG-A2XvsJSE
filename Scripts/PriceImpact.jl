## Author: Patrick Chang
# Script file to plot the price impact of 10 equities
# from JSE and A2X over the same period of 2019-01-02 to 2019-07-15.

#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots
using Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

# JSE tickers
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]

# A2X tickers
A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]

# Load the price impact data from each exchange

A2X_PriceImpact = load("Real Data/A2X/PriceImpact/A2X_PriceImpact.jld")
A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Real Data/JSE/PriceImpact/JSE_PriceImpact.jld")
JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]

## Function to construct the plotting information for impact
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

## Function to plot the price impact
function PlotImpact(data, ticker::Vector, side::Symbol; low = -3, up = 1)
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
        push!(Impact, i => getPriceImpact(data[i], low = low, up = up))
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
        xlabel!(L"\textrm{Buyer-Initiated: } \omega^*")
        ylabel!(L"\Delta p^*")
    elseif side == :sell
        xlabel!(L"\textrm{Seller-Initiated: } \omega^*")
        ylabel!(L"\Delta p^*")
    end
end


## Obtain results

PlotImpact(JSE_PriceImpact, JSE_tickers, :buy, low = -3, up = 1)
# savefig("Plots/JSEImpactBuy.svg")

PlotImpact(JSE_PriceImpact, JSE_tickers, :sell, low = -3, up = 1)
# savefig("Plots/JSEImpactSell.svg")

PlotImpact(A2X_PriceImpact, A2X_tickers, :buy, low = -3, up = 1)
# savefig("Plots/A2XImpactBuy.svg")

PlotImpact(A2X_PriceImpact, A2X_tickers, :sell, low = -3, up = 1)
# savefig("Plots/A2XImpactSell.svg")
