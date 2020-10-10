## Author: Patrick Chang and Ivan Jerivich
# Script file to investigate the accuracy of rules to infer
# the trade classification.
# Here we only perform this for A2X data, as we have unique access
# to the true trade classification in this case.

## Preamble

using CSV, DataFrames, JLD, Dates, ProgressMeter, Plots, Statistics, LaTeXStrings

cd("/Users/patrickchang1/PCIJAPTG-A2XvsJSE")

A2X = CSV.read("Real Data/A2X/Cleaned/A2X_Cleaned_NPN.csv")
