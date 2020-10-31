### Title: Price impact
### Authors: Patrick Chang and Ivan Jericevich
### Function: Visualise the price impact curves of JSE and A2X equities over the same period
### Structure:
# 1. Preliminaries
# 2. Compute price impact curves
# 3. Visualization
#---------------------------------------------------------------------------


### 1. Preamble
using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------


### 2. Compute price impact curves
function getPriceImpact(data::DataFrame; low = -3, up = 1) # Function to construct the plotting information for impact
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
#---------------------------------------------------------------------------


### 3. Visualization
PlotImpact(JSE_PriceImpact, JSE_tickers, :buy, low = -3, up = 1); savefig("Figures/JSEImpactBuy.pdf")
PlotImpact(JSE_PriceImpact, JSE_tickers, :sell, low = -3, up = 1); savefig("Figures/JSEImpactSell.pdf")
PlotImpact(A2X_PriceImpact, A2X_tickers, :buy, low = -3, up = 1); savefig("Figures/A2XImpactBuy.pdf")
PlotImpact(A2X_PriceImpact, A2X_tickers, :sell, low = -3, up = 1); savefig("Figures/A2XImpactSell.pdf")
#---------------------------------------------------------------------------
