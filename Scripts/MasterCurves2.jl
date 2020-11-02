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
using CSV, DataTables, DataFrames, JLD, Dates, ProgressMeter, Plots, Optim, Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
JSE_tickers = ["ABG", "AGL", "BTI", "FSR", "NED", "NPN", "SBK", "SHP", "SLM", "SOL"]; A2X_tickers = ["APN", "ARI", "AVI", "CML", "GRT", "MRP", "NPN", "SBK", "SLM", "SNT"]
A2X_PriceImpact = load("Test Data/A2X/Price Impact/A2X_PriceImpact.jld"); A2X_PriceImpact = A2X_PriceImpact["A2X_PriceImpact"]
JSE_PriceImpact = load("Test Data/JSE/Price Impact/JSE_PriceImpact.jld"); JSE_PriceImpact = JSE_PriceImpact["JSE_PriceImpact"]
#---------------------------------------------------------------------------

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
function Collapse(side::Symbol, paramJSE::Vector, paramA2X::Vector, dataJSE, dataA2X, tickerJSE = JSE_tickers, tickerA2X = A2X_tickers, low = -1, up = 1) # Plot the Master impact from each exchange together
    # Extract liquidity
    ADVJSE = dataJSE[3]; ADVA2X = dataA2X[3]
    # Extract appropriate side
    if side == :buy
        # Buyer-Initiated
        dataJSE = dataJSE[1]; dataA2X = dataA2X[1]
    elseif side == :sell
        # Seller-Initiated
        dataJSE = dataJSE[2]; dataA2X = dataA2X[2]
    end
    # Get Scaled impact from JSE
    JSEΔp = fill(NaN, 20, 10); JSEω = fill(NaN, 20, 10)
    for i in 1:length(tickerJSE)
        temp = getMasterImpact(dataJSE[tickerJSE[i]], ADVJSE[tickerJSE[i]], paramJSE, low = low, up = up)
        inds = temp[3]
        JSEΔp[inds, i] = temp[2]; JSEω[inds, i] = temp[1]
    end
    # Get Scaled impact from A2X
    A2XΔp = fill(NaN, 20, 10); A2Xω = fill(NaN, 20, 10)
    scalingParameter = side == :sell ? 1 : 1 # 9000 : 1
    for i in 1:length(tickerA2X)
        temp = getMasterImpact(dataA2X[tickerA2X[i]], ADVA2X[tickerA2X[i]], paramA2X, low = low, up = up)
        inds = temp[3]
        z = temp[2] ./ scalingParameter
        A2XΔp[inds, i] = z
        A2Xω[inds, i] = temp[1]
    end
    # Get the mean and std
    collapsedJSEω = mean(JSEω,dims=2)
    collapsedJSEΔp = mean(JSEΔp,dims=2)
    varJSEΔp = std(JSEΔp,dims=2)
    collapsedA2Xω = fill(0.0, 20, 1)
    collapsedA2XΔp = fill(0.0, 20, 1)
    varA2XΔp = fill(0.0, 20, 1)
    for i in 1:20
        collapsedA2Xω[i] = mean(filter(!isnan, A2Xω[i,:]))
        tempΔp = filter(!isnan, A2XΔp[i,:])
        collapsedA2XΔp[i] = mean(tempΔp)
        varA2XΔp[i] = std(tempΔp)
    end
    return collapsedA2Xω, collapsedA2XΔp, collapsedJSEω, collapsedJSEΔp
end
function Bootstrap(M::Int64, JSEData::Tuple{Dict{Any,Any},Dict{Any,Any},Dict{Any,Any}}, A2XData::Tuple{Dict{Any,Any},Dict{Any,Any},Dict{Any,Any}}, initiated::Symbol)
    # Get parameters and tickers
    if initiated == :Buy
        parametersJSE = load("Test Data/JSE/Price Impact/JSEParams.jld")["JSEBuyParam"]; parametersA2X = load("Test Data/A2X/Price Impact/A2XParams.jld")["A2XBuyParam"] # Load scaling parameters
        tickersJSE = keys(JSEData[1]); tickersA2X = keys(A2XData[1]) # Extract tickers
    else
        parametersJSE = load("Test Data/JSE/Price Impact/JSEParams.jld")["JSESellParam"]; parametersA2X = load("Test Data/A2X/Price Impact/A2XParams.jld")["A2XSellParam"] # Load scaling parameters
        tickersJSE = keys(JSEData[2]); tickersA2X = keys(A2XData[2]) # Extract tickers
    end
    # Initialise storage for bootstrap master curves
    JSEω = fill(0.0, (M, 20)); JSEΔp = fill(0.0, (M, 20)); A2Xω = fill(0.0, (M, 20)); A2XΔp = fill(0.0, (M, 20))
    @showprogress "Computing:" for m in 1:M # Number of bootstrap samples
        if initiated == :Buy
            # Extract relevant side
            JSE = JSEData[1]; A2X = A2XData[1]
            # Initialise dictionary of bootstrap samples with keys corresponding to tickers
            bootstrapSampleJSE = Dict(); bootstrapSampleA2X = Dict()
            # Iterate through tickers in each exchange and create bootstrap samples from price impact data
            for (JSETicker, A2XTicker) in zip(tickersJSE, tickersA2X)
                # Remove NaNs in price impact data
                JSETemp = filter(x -> !isnan(x.Impact), JSE[JSETicker]); A2XTemp = filter(x -> !isnan(x.Impact), A2X[A2XTicker])
                # Sample with replacement
                bootstrapIndecesJSE = sample(1:nrow(JSETemp), nrow(JSETemp), replace = true); bootstrapIndecesA2X = sample(1:nrow(A2XTemp), nrow(A2XTemp), replace = true) # Sample with replacement to get bootstrap sample
                bootstrapSampleJSE[JSETicker] = JSETemp[bootstrapIndecesJSE, :]; bootstrapSampleA2X[A2XTicker] = A2XTemp[bootstrapIndecesA2X, :]
            end
            # Craete new bootstrapped price impact data
            JSEPriceImpact = (bootstrapSampleJSE, JSEData[2], JSEData[3]); A2XPriceImpact = (bootstrapSampleA2X, A2XData[2], A2XData[3])
            # Construct collapsed master curves from bootstrap sample
            collapsed = Collapse(:buy, parametersJSE, parametersA2X, JSEPriceImpact, A2XPriceImpact)
        else
            # Extract relevant side
            JSE = JSEData[2]; A2X = A2XData[2]
            # Initialise dictionary of bootstrap samples with keys corresponding to tickers
            bootstrapSampleJSE = Dict(); bootstrapSampleA2X = Dict()
            # Iterate through tickers in each exchange and create bootstrap samples from price impact data
            for (JSETicker, A2XTicker) in zip(tickersJSE, tickersA2X)
                # Remove NaNs in price impact data
                JSETemp = filter(x -> !isnan(x.Impact), JSE[JSETicker]); A2XTemp = filter(x -> !isnan(x.Impact), A2X[A2XTicker])
                # Sample with replacement
                bootstrapIndecesJSE = sample(1:nrow(JSETemp), nrow(JSETemp), replace = true); bootstrapIndecesA2X = sample(1:nrow(A2XTemp), nrow(A2XTemp), replace = true) # Sample with replacement to get bootstrap sample
                bootstrapSampleJSE[JSETicker] = JSETemp[bootstrapIndecesJSE, :]; bootstrapSampleA2X[A2XTicker] = A2XTemp[bootstrapIndecesA2X, :]
            end
            # Craete new bootstrapped price impact data
            JSEPriceImpact = (JSEData[1], bootstrapSampleJSE, JSEData[3]); A2XPriceImpact = (A2XData[1], bootstrapSampleA2X, A2XData[3])
            # Construct collapsed master curves from bootstrap sample
            collapsed = Collapse(:sell, parametersJSE, parametersA2X, JSEPriceImpact, A2XPriceImpact)
        end
        # Store
        A2Xω[m, :] = collapsed[1]; A2XΔp[m, :] = collapsed[2]; JSEω[m, :] = collapsed[3]; JSEΔp[m, :] = collapsed[4]
    end
    return A2Xω, A2XΔp, JSEω, JSEΔp
end
bootstrapBuy = Bootstrap(100, JSE_PriceImpact, A2X_PriceImpact, :Buy); bootstrapSell = Bootstrap(100, JSE_PriceImpact, A2X_PriceImpact, :Sell)
function PlotJoint(p, M::Int64, initiated::Symbol; JSEPriceImpact = JSE_PriceImpact, A2XPriceImpact = A2X_PriceImpact)
    bootstrap = Bootstrap(M, JSEPriceImpact, A2XPriceImpact, initiated)
    upperA2X = Vector{Float64}(); upperJSE = Vector{Float64}(); lowerA2X = Vector{Float64}(); lowerJSE = Vector{Float64}(); averageA2XΔp = Vector{Float64}();  averageA2Xω = Vector{Float64}(); averageJSEΔp = Vector{Float64}(); averageJSEω = Vector{Float64}()
    for j in 1:size(bootstrap[2])[2]
        # Extract non NaN indeces in the normalised volumes for the x-axis of mean curves and error bars
        indecesA2Xω = findall(x -> !isnan(x), bootstrap[1][:, j]); indecesJSEω = findall(x -> !isnan(x), bootstrap[3][:, j])
        # Extract non NaN indeces in the impacts to get the mean impacts and error bars
        indecesA2XΔp = findall(x -> !isnan(x), bootstrap[2][:, j]); indecesJSEΔp = findall(x -> !isnan(x), bootstrap[4][:, j])
        if !isempty(indecesA2Xω) && !isempty(indecesA2XΔp)
            push!(averageA2Xω, mean(bootstrap[1][indecesA2Xω, j])); push!(averageA2XΔp, mean(bootstrap[2][indecesA2XΔp, j]))
            push!(upperA2X, quantile(bootstrap[2][indecesA2XΔp, j], 0.95)); push!(lowerA2X, quantile(bootstrap[2][indecesA2XΔp, j], 0.05))
            # push!(upperA2X, averageA2XΔp[end] + std(bootstrap[2][indecesA2XΔp, j])); push!(lowerA2X, averageA2XΔp[end] - std(bootstrap[2][indecesA2XΔp, j]))
        end
        if !isempty(indecesJSEω) && !isempty(indecesJSEΔp)
            push!(averageJSEω, mean(bootstrap[3][indecesJSEω, j])); push!(averageJSEΔp, mean(bootstrap[4][indecesJSEΔp, j]))
            push!(upperJSE, quantile(bootstrap[4][indecesJSEΔp, j], 0.95)); push!(lowerJSE, quantile(bootstrap[4][indecesJSEΔp, j], 0.05))
            # push!(upperJSE, averageJSEΔp[end] + std(bootstrap[4][indecesJSEΔp, j])); push!(lowerJSE, averageJSEΔp[end] - std(bootstrap[4][indecesJSEΔp, j]))
        end
    end
    # Plot results
    if initiated == :Buy
        plot!(p, averageA2Xω, averageA2XΔp, linecolor = :blue, fillalpha = .3, fillcolor = :blue, seriestype = :line, ribbon = (lowerA2X, upperA2X), scale = :log10, label = "")
        plot!(p, averageA2Xω, averageA2XΔp, markercolor = :blue, markerstrokecolor = :blue, seriestype = :scatter, scale = :log10, label = "extrm{A2X - Buyer Initiated}")
        plot!(p, averageJSEω, averageJSEΔp, linecolor = :red, fillalpha = .3, fillcolor = :red, seriestype = :line, ribbon = (lowerJSE, upperJSE), scale = :log10, label = "")
        plot!(p, averageJSEω, averageJSEΔp, markercolor = :red, markerstrokecolor = :red, seriestype = :scatter, scale = :log10, label = "extrm{JSE - Buyer Initiated}")
    else
        plot!(p, averageA2Xω, averageA2XΔp, linecolor = :blue, fillalpha = .3, fillcolor = :blue, seriestype = :line, ribbon = (lowerA2X, upperA2X), scale = :log10, label = "")
        plot!(p, averageA2Xω, averageA2XΔp, markercolor = :blue, markerstrokecolor = :blue, markershape = :utriangle, seriestype = :scatter, scale = :log10, label = "extrm{A2X - Seller Initiated}")
        plot!(p, averageJSEω, averageJSEΔp, linecolor = :red, fillalpha = .3, fillcolor = :red, seriestype = :line, ribbon = (lowerJSE, upperJSE), scale = :log10, label = "")
        plot!(p, averageJSEω, averageJSEΔp, markercolor = :red, markerstrokecolor = :red, markershape = :utriangle, seriestype = :scatter, scale = :log10, label = "extrm{JSE - Seller Initiated}")
    end
end



p = plot(xlabel = L"\overline{\omega^* / C^{\delta}}", ylabel = L"\overline{\Delta p^* C^{\gamma}} \kappa", legend = :outertopright, legendtitle = L"\textrm{Ticker}", size = (700, 400), dpi = 300)
PlotJoint(p, 100, :Sell); PlotJoint(p, 100, :Buy)
savefig("z.png")
