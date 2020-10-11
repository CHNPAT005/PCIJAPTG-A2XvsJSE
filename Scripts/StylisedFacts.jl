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
    filter!(x-> x >= Date("2019-01-01"), dates_unique)
    # Create master dataframe to store all useful information
    Rets = Float64[]
    # Loop through each day and piece the returns from there
    @showprogress "Computing..." for k in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[k]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Pull out the item to make returns with
        prices = tempdata[item]
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
function PLqqplot_Tail(obs, side) # Function to plot the QQ plot using the 95% quantile as x_min
    # Order and sort observations
    nobs = length(obs)
    obs = side == "Left" ? sort(-obs) : sort(obs)
    # Find 95% quantile and filter out data from there
    q95 = Int(round(0.95*nobs))
    obs = obs[q95:end]
    # Get appropriate Xmin and α
    xmin = obs[1]
    α = PLalpha(obs, xmin)
    # Get quantiles
    nobs = length(obs)
    quantiles⁰ = [PLquantile(i/nobs, α, xmin) for i in 1:nobs]
    # Density + QQ plot
    density(getTail(obs, "Left"), label = "", color = :blue)
    xlabel!(L"\textrm{Price fluctuations}")
    ylabel!(L"\textrm{Density}")
    # Plot the QQ plot
    plot!([quantiles⁰ quantiles⁰], [obs quantiles⁰], seriestype = [:scatter :line],
    inset = (1, bbox(0.35, 0.1, 0.65, 0.65, :top)), subplot = 2, legend = :none, xlabel = L"\textrm{Theoretical Quantiles}", ylabel = L"\textrm{Sample Quantiles}", title = L"\textrm{Power Law QQ-plot}", guidefontsize = 8, linecolor = :black, markercolor = :blue, markerstrokecolor = :blue)
end
#---------------------------------------------------------------------------


# 4. Generate the stylized facts
for interval in [1, 10, 20, 30]
    ## JSE
    NPNJSEMicroBars = CSV.read(string("Real Data/JSE/Bar/NPNMicroPriceBars", interval, "min.csv"))
    NPNJSEMicroReturns = getReturns(NPNJSEMicroBars, :Close)
    # Return auto-correlations
    ACFJSE = plot(1:100, autocor(convert.(Float64, NPNJSEMicroReturns), 1:100), label = "", seriestype = :sticks, xlabel = L"\textrm{Lag}", ylabel = L"\textrm{ACF}", color = :black)
    hline!(ACFJSE, [quantile(Normal(), (1 + 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :blue, label = "")
    hline!(ACFJSE, [quantile(Normal(), (1 - 0.95) / 2) / sqrt(length(NPNJSEMicroReturns))], color = :blue, label = "")
    savefig(ACFJSE, string("Plots/Returns ACF ", interval, "min (JSE).pdf"))
    # Density plots of right and left tails of returns
    rightTailDensityJSE = density(getTail(NPNJSEMicroReturns, "Right"), label = "")
    xlabel!(rightTailDensityJSE, L"\textrm{Price fluctuations}")
    ylabel!(rightTailDensityJSE, L"\textrm{Density}")
    savefig(rightTailDensityJSE, string("Plots/Right Tail Returns Density ", interval, "min (JSE).pdf"))
    leftTailDensityJSE = density(getTail(NPNJSEMicroReturns, "Left"), label = "")
    xlabel!(leftTailDensityJSE, L"\textrm{Price fluctuations}")
    ylabel!(leftTailDensityJSE, L"\textrm{Density}")
    savefig(leftTailDensityJSE, string("Plots/Left Tail Returns Density ", interval, "min (JSE).pdf"))
    # QQ-plots of right and left tails of returns relative to a fitted power-law distribution
    rightTailQQJSE = PLqqplot_Tail(NPNJSEMicroReturns, "Right")
    savefig(rightTailQQJSE, string("Plots/Right Tail Returns QQ-plot ", interval, "min (JSE).pdf"))
    leftTailQQJSE = PLqqplot_Tail(NPNJSEMicroReturns, "Left")
    savefig(leftTailQQJSE, string("Plots/Left Tail Returns QQ-plot ", interval, "min (JSE).pdf"))
    ## A2X
end
#---------------------------------------------------------------------------
