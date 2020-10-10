## Author: Patrick Chang and Ivan Jerivich
# Script file to look at bar data stylised facts between A2X and JSE.
# 1) Visualise micro-price bars
# 2) Micro-price bar returns for the 95% tails
# 3) Micro-price bar return auto-correlation
# Note that we investigate three time-scales: 1min, 10min and 20min

## Preamble

using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

A2X = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")
JSE = CSV.read("Real Data/JSE/Cleaned/JSECleanedTAQNPN.csv")
