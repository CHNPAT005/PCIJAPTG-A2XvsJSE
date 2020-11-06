# PCIJAPTG-A2XvsJSE

<img src="Figures/A2XLogo.png" width="48"> <img src="Figures/JSELogo.png" width="48">

## Authors
* Ivan Jericevich
* Patrick Chang
* Aveshen Pillay
* Tim Gebbie

## Description
This repository details the implementation of an exploratory data analysis of the South African market microstructure. Specifically, it considers shared listings on two South African equity exchanges: the Johannesburg Stock Exchange and the A2X exchange. The related paper explores an empircal comparisons of markets trading similar shares, in a similar regulatory and economic environment, but with vastly different liquidity, cost and business models. We compare the distributions and auto-correlations of returns on different time scales, we compare price impact and master curves, and we compare the cost of trading on each exchange.

|A2X                     |                 |JSE                         |                 |
|------------------------|-----------------|----------------------------|-----------------|
|**Security name**       |**Security code**|**Security name**           |**Security code**|
|Aspen Pharmacare        |APN              |Absa Group Ltd              |ABG              |
|African Rainbow Min Ltd |ARI              |Anglo American Plc          |AGL              |
|AVI Ltd                 |AVI              |British American Tobacco Plc|BTI              |
|Coronation Fund Managers|CML              |FirstRand Limited           |FSR              |
|Growthpoint Prop Ltd    |GRT              |Nedbank Group Ltd           |NED              |
|Mr Price Group Limited  |MRP              |Naspers Ltd                 |NPN              |
|Naspers Ltd             |NPN              |Standard Bank Group Ltd     |SBK              |
|Standard Bank Group Ltd |SBK              |Shoprite Holdings Ltd       |SHP              |
|Sanlam Limited          |SLM              |Sanlam Limited              |SLM              |
|Santam Limited          |SNT              |Sasol Ltd                   |SOL              |


## Prerequisites
* [Julia](https://julialang.org) programming langauge
* A text editor or IDE such as [Atom](https://flight-manual.atom.io/getting-started/sections/installing-atom/), [VS Code](https://code.visualstudio.com/download) or [Jupyter](https://jupyter.org/install)
* The test data to be used for the replication of results and algorithm validation can be found at [ZivaHub](https://figshare.com/articles/dataset/_/13187591)

## Usage
Clone the repository
```sh
git clone https://github.com/CHNPAT005/PCIJAPTG-A2XvsJSE.git
```
and set the working directory in all `.jl` files.
```julia
cd("C:/Users/.../PCIJAPTG-A2XvsJSE")
```
Packages can be installed using
```julia
Pkg.add("...")
```
The datasets used for the paper is for the period between 2019-01-02 to 2019-07-15. The data is proprietary and is not included here. For the purpose of algorithm verification we include the last week of data from both exchanges (2019-07-8 - 2019-07-12). The data for this can be downloaded [here](https://figshare.com/articles/dataset/_/13187591). There should be 5 A2X files and 10 JSE files. These should be placed in `Test Data/A2X/Raw/` and `Test Data/JSE/Raw/` respectively. The Computed data folder is in place to hold computed results (`.jld` files) that apply to both exchanges. All other data specific to each exchange will be saved to the Test Data folder (upon implementation). All figures contained within the Figures folder correspond to those in the paper and were obtained from the implementation on the full data sets. Data cleaning scripts should be run first before any of the Julia files in the Scripts folder can be run.
