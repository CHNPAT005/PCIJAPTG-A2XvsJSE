## Author: Patrick Chang and Ivan Jerivich
# Script file to investigate:
# 1) Order flow auto-correlation on A2X and JSE
# 2) Micro-price return auto-correlation on A2X and JSE

## Preamble

using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

A2X = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")
JSE = CSV.read("Real Data/JSE/Cleaned/JSECleanedTAQNPN.csv")
