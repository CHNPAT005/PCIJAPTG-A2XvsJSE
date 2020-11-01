# PCIJAPTG-A2XvsJSE

<img src="Figures/A2XLogo.png" width="48"> <img src="Figures/JSELogo.png" width="48">

## Authors
* Ivan Jericevich
* Patrick Chang
* Aveshen Pillay
* Tim Gebbie

## Introduction
This repository details the implementation of an exploratory data analysis of the South African market microstructure. Specifically, it considers shared listings on two South African equity exchanges: the Johannesburg Stock Exchange and the A2X exchange.

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
* [Julia](http://movielens.org) programming langauge
* A text editor or IDE such as [Atom](https://flight-manual.atom.io/getting-started/sections/installing-atom/), [VS Code](https://code.visualstudio.com/download) or [Jupyter](https://jupyter.org/install)
* The test data sets accompanied herewith which can also be found at <http://linktozivahub/datasets/>

## Usage
Clone the repository
```sh
git clone https://github.com/CHNPAT005/PCIJAPTG-A2XvsJSE.git
```
The datasets used to produce all figures from both exchanges is for the period 2019-01-01 - 2019-07-15. This data is proprietary and is thus not included here. For the purpose of replication and proof of concept, we include the last week of data from both exchanges (2019-07-8 - 2019-07-12) in the [Test Data](..tree/release/Test Data) folder.