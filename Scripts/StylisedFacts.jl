### Title: Stylized Facts
### Authors: Patrick Chang and Ivan Jericevich
### Function: Obtain the figures which exhibit the intra-day stylized facts for Naspers in both the JSE and A2X market
### Structure:
# 1. Preliminaries
# 2. Functions for generating the upper and lower 95-percentiles of micro-price returns
# 3. Functions for fitting a power-law distribution and generating a QQ-plot of the upper and lower tails of micro-price returns
# 4. Generate the stylized facts
#   - Auto-correlations of micro-price returns
#   - Densities of the upper and lower 95-percentiles of micro-price returns
#   - QQ-plots of right and left tails of returns relative to a fitted power-law distribution
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, DataFrames, Dates, ProgressMeter, Plots, StatsPlots, Statistics, LaTeXStrings, TimeSeries, Distributions, StatsBase
cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")
#---------------------------------------------------------------------------


### 2. Functions for generating the upper and lower 95-percentiles of micro-price returns
function getReturns(data::DataFrame, item::Symbol)
    # Extract the dates
    dates = Date.(data[:,1])
    # Get unique dates and start from 2019-01-01 onwards
    dates_unique = unique(dates)
    filter!(x-> Date("2019-01-01") <= x <= Date("2019-07-15"), dates_unique) ################################ Fixed error
    # Create master dataframe to store all useful information
    Rets = Float64[]
    # Loop through each day and piece the returns from there
    @showprogress "Computing..." for k in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[k]
        tempdata = data[findall(x -> x == tempday, dates), :]
        filter!(x -> Time(x.TimeStamp) > Time(9, 01, 00), tempdata)
        # Pull out the item to make returns with
        prices = tempdata[:, item]
        # Returns are price fluctuations
        ret = diff(log.(filter(!isnan, prices)))
        # Append to vector of returns
        Rets = [Rets; ret]
    end
    return Rets
end
function getTail(data::Vector, side::String)
    # Get length of data
    nobs = length(data)
    # Sort the data in reverse
    data = side == "Left" ? sort(-data) : sort(data)
    # Find indicies
    q95 = Int(round(0.95*nobs))
    # 95% and above of negative of sorted values
    data = data[q95:end]
    return data
end
#---------------------------------------------------------------------------


### 3. Functions for fitting a power-law distribution and generating a QQ-plot of the upper and lower tails of micro-price returns
function PLalpha(data, xmin) # Function to compute the MLE for α in the power law distribution. Input parameters are full data + x_min
    data = filter(x -> x >= xmin, data)
    α = 1 + length(data) / sum(log.(data ./ xmin))
    return α
end
function PLquantile(p, α, xmin) # Function to get the quantile for Power Law
    # From the CDF, obtain the quantile
    return (1 - p)^(-1 / (α - 1)) * xmin
end
#---------------------------------------------------------------------------
################ Joint Plots ######################

function PLqqplot_Tail(p, obs, side, market) # Function to plot the QQ plot using the 95% quantile as x_min
    # Order and sort observations
    obs = getTail(obs, side) ############################# Fixed error
    # Get appropriate Xmin and α
    xmin = obs[1]
    α = PLalpha(obs, xmin)
    # Get quantiles
    nobs = length(obs)
    quantiles⁰ = [PLquantile(i/nobs, α, xmin) for i in 1:nobs]
    # Density + QQ plot
    if market == "A2X"
        col = :blue
        pos = (1, bbox(0.6, 0.03, 0.34, 0.34, :top))
        sub = 3
        density!(p, obs, label = "A2X", color = col, legend = :topleft)
    else
        col = :red
        pos = (1, bbox(0.6, 0.5, 0.34, 0.34, :top))
        sub = 2
        density!(p, obs, label = "JSE", color = col, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}", legend = :topleft)
    end

    # Plot the QQ plot
    plot!(p, [quantiles⁰ quantiles⁰], [obs quantiles⁰], seriestype = [:scatter :line], inset = pos, subplot = sub, legend = :none, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{%$market - Power Law QQ-plot}", titlefontsize = 7, guidefontsize = 7, tickfont = 5, linecolor = :black, markercolor = col, markerstrokecolor = col)
end

for interval in [1, 10, 20, 30]
    NPNJSEMicroBars = CSV.read(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNJSEMicroReturns = getReturns(NPNJSEMicroBars, :Close)
    NPNA2XMicroBars = CSV.read(string("Real Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNA2XMicroReturns = getReturns(NPNA2XMicroBars, :Close)
    # Return auto-correlations
    ACFJSE = plot(1:100, autocor(convert.(Float64, NPNJSEMicroReturns), 1:100), label = "", seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", color = :black, dpi = 300)
    hline!(ACFJSE, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :red, label = "")
    hline!(ACFJSE, [quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :red, label = "")
    savefig(ACFJSE, string("Plots/ReturnsACFJSE", interval, "min.svg"))
    ACFA2X = plot(1:100, autocor(convert.(Float64, NPNA2XMicroReturns), 1:100), label = "", seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", color = :black, dpi = 300)
    hline!(ACFA2X, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(NPNA2XMicroReturns))], color = :blue, label = "")
    hline!(ACFA2X, [quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(NPNA2XMicroReturns))], color = :blue, label = "")
    savefig(ACFA2X, string("Plots/ReturnsACFA2X", interval, "min.svg"))
    # Joint Density and QQ-plots of right and left tails of returns relative to a fitted power-law distribution
    rightTailQQ = plot(dpi = 300)
    PLqqplot_Tail(rightTailQQ, NPNJSEMicroReturns, "Right", "JSE")
    PLqqplot_Tail(rightTailQQ, NPNA2XMicroReturns, "Right", "A2X")
    savefig(rightTailQQ, string("Plots/RightTailReturns", interval, "min.svg"))
    leftTailQQ = plot(dpi = 300)
    PLqqplot_Tail(leftTailQQ, NPNJSEMicroReturns, "Left", "JSE")
    PLqqplot_Tail(leftTailQQ, NPNA2XMicroReturns, "Left", "A2X")
    savefig(leftTailQQ, string("Plots/LeftTailReturns", interval, "min.svg"))
    # Joint density plot of returns
    densityReturns = density(NPNJSEMicroReturns, label = "JSE", color = :red, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}", legend = :bottomright)
    density!(densityReturns, NPNA2XMicroReturns, label = "A2X", color = :blue, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}")
    qqplot!(densityReturns, Normal, NPNJSEMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{JSE - Normal QQ-plot}", markersize = 3, markercolor = :red, markerstrokecolor = :red, linecolor = :black, inset = (1, bbox(0.62, 0.1, 0.33, 0.33, :top)), subplot = 2, legend = :none, titlefontsize = 7, guidefontsize = 7, tickfont = 5)
    qqplot!(densityReturns, Normal, NPNA2XMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{A2X - Normal QQ-plot}", markersize = 3, markercolor = :blue, markerstrokecolor = :blue, linecolor = :black, inset = (1, bbox(0.1, 0.1, 0.33, 0.33, :top)), subplot = 3, legend = :none, titlefontsize = 7, guidefontsize = 7, tickfont = 5)
    savefig(densityReturns, string("Plots/ReturnsDistribution ", interval, "min.pdf"))
end


# Tick-by-tick return distributions
NPNJSE = CSV.read("Real Data/JSE/JSECleanedTAQ" * "NPN" * ".csv")
NPNA2X = CSV.read("Real Data/A2X/A2X_Cleaned_" * "NPN" * ".csv")[:, 2:end]
DataFrames.rename!(NPNA2X, (:Date => :TimeStamp))
NPNJSEMicroReturns = getReturns(NPNJSE, :MicroPrice); NPNA2XMicroReturns = getReturns(NPNA2X, :MicroPrice)

rightTailQQ = plot(dpi = 300)
PLqqplot_Tail(rightTailQQ, NPNJSEMicroReturns, "Right", "JSE")
PLqqplot_Tail(rightTailQQ, NPNA2XMicroReturns, "Right", "A2X")
savefig(rightTailQQ, string("Plots/RightTailTickReturns.svg"))
leftTailQQ = plot(dpi = 300)
PLqqplot_Tail(leftTailQQ, NPNJSEMicroReturns, "Left", "JSE")
PLqqplot_Tail(leftTailQQ, NPNA2XMicroReturns, "Left", "A2X")
savefig(leftTailQQ, string("Plots/LeftTailTickReturns.svg"))
# Joint density plot of returns
densityReturns = density(NPNJSEMicroReturns, label = "JSE", color = :red, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}", legend = :bottomright)
density!(densityReturns, NPNA2XMicroReturns, label = "A2X", color = :blue, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}")
qqplot!(densityReturns, Normal, NPNJSEMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{JSE - Normal QQ-plot}", markersize = 3, markercolor = :red, markerstrokecolor = :red, linecolor = :black, inset = (1, bbox(0.62, 0.1, 0.33, 0.33, :top)), subplot = 2, legend = :none, titlefontsize = 7, guidefontsize = 7, tickfont = 5)
qqplot!(densityReturns, Normal, NPNA2XMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{A2X - Normal QQ-plot}", markersize = 3, markercolor = :blue, markerstrokecolor = :blue, linecolor = :black, inset = (1, bbox(0.1, 0.1, 0.33, 0.33, :top)), subplot = 3, legend = :none, titlefontsize = 7, guidefontsize = 7, tickfont = 5)
savefig(densityReturns, string("Plots/TickReturnsDistribution.pdf"))










#---------------------------------------------------------------------------

################ Disjoint Plots ######################

function PLqqplot_Tail(p, obs, side, market) # Function to plot the QQ plot using the 95% quantile as x_min
    # Order and sort observations
    obs = getTail(obs, side)
    # Get appropriate Xmin and α
    xmin = obs[1]
    α = PLalpha(obs, xmin)
    # Get quantiles
    nobs = length(obs)
    quantiles⁰ = [PLquantile(i/nobs, α, xmin) for i in 1:nobs]
    # Density + QQ plot
    col = market == "A2X" ? :blue : :red
    pos = market == "A2X" ? :blue : :red
    density!(p, obs, label = "", color = col, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}")
    # Plot the QQ plot
    plot!(p, [quantiles⁰ quantiles⁰], [obs quantiles⁰], seriestype = [:scatter :line],
    inset = (1, bbox(0.6, 0.1, 0.35, 0.35, :top)), subplot = 2, legend = :none, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{Power Law QQ-plot}", guidefontsize = 8, linecolor = :black, markercolor = col, markerstrokecolor = col)
    return p
end


# 4. Generate the stylized facts
for interval in [1, 10, 20, 30]
    ## JSE
    NPNJSEMicroBars = CSV.read(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNJSEMicroReturns = getReturns(NPNJSEMicroBars, :Close)
    # Return auto-correlations
    ACFJSE = plot(1:100, autocor(convert.(Float64, NPNJSEMicroReturns), 1:100), label = "", seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", color = :black)
    hline!(ACFJSE, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :red, label = "")
    hline!(ACFJSE, [quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :red, label = "")
    savefig(ACFJSE, string("Plots/Returns ACF ", interval, "min (JSE).pdf"))
    # Density and QQ-plots of right and left tails of returns relative to a fitted power-law distribution
    rightTailQQJSE = PLqqplot_Tail(NPNJSEMicroReturns, "Right", "JSE")
    savefig(rightTailQQJSE, string("Plots/Right Tail Returns QQ-plot ", interval, "min (JSE).pdf"))
    leftTailQQJSE = PLqqplot_Tail(NPNJSEMicroReturns, "Left", "JSE")
    savefig(leftTailQQJSE, string("Plots/Left Tail Returns QQ-plot ", interval, "min (JSE).pdf"))
    # Density plots of returns
    densityJSE = density(NPNJSEMicroReturns, label = "", color = :red, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}")
    qqplot!(densityJSE, Normal, NPNJSEMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", markersize = 3, markercolor = :red, markerstrokecolor = :red, linecolor = :black, inset = (1, bbox(0.6, 0.1, 0.35, 0.35, :top)), subplot = 2, legend = :none,  guidefontsize = 8)
    savefig(densityJSE, string("Plots/Returns Distribution ", interval, "min (JSE).pdf"))
    ## A2X
    NPNA2XMicroBars = CSV.read(string("Real Data/A2X/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNA2XMicroReturns = getReturns(NPNA2XMicroBars, :Close)
    # Return auto-correlations
    ACFA2X = plot(1:100, autocor(convert.(Float64, NPNA2XMicroReturns), 1:100), label = "", seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", color = :black)
    hline!(ACFA2X, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(NPNA2XMicroReturns))], color = :blue, label = "")
    hline!(ACFA2X, [quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(NPNA2XMicroReturns))], color = :blue, label = "")
    savefig(ACFA2X, string("Plots/Returns ACF ", interval, "min (A2X).pdf"))
    # Density and QQ-plots of right and left tails of returns relative to a fitted power-law distribution
    rightTailQQA2X = PLqqplot_Tail(NPNA2XMicroReturns, "Right", "A2X")
    savefig(rightTailQQA2X, string("Plots/Right Tail Returns QQ-plot ", interval, "min (A2X).pdf"))
    leftTailQQA2X = PLqqplot_Tail(NPNA2XMicroReturns, "Left", "A2X")
    savefig(leftTailQQA2X, string("Plots/Left Tail Returns QQ-plot ", interval, "min (A2X).pdf"))
    # Density plots of returns
    densityA2X = density(NPNA2XMicroReturns, label = "", color = :blue, xlabel = L"\textrm{Price fluctuations}", ylabel = L"\textrm{Density}")
    qqplot!(densityA2X, Normal, NPNA2XMicroReturns, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", markersize = 3, markercolor = :blue, markerstrokecolor = :blue, linecolor = :black, inset = (1, bbox(0.6, 0.1, 0.35, 0.35, :top)), subplot = 2, legend = :none,  guidefontsize = 8)
    savefig(densityA2X, string("Plots/Returns Distribution ", interval, "min (A2X).pdf"))
end
#---------------------------------------------------------------------------
