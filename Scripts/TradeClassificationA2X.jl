### Title: Trade Classification - A2X
### Authors: Patrick Chang and Ivan Jericevich
### Function: Investigate the accuracy of rules to infer the trade classification
### Structure:
# 1. Preliminaries
# 2. Classification rules
# 3. Implement classification rules
#---------------------------------------------------------------------------


### 1. Preliminaries
using CSV, CodecBzip2, DataFrames, ProgressMeter, Dates, JLD, Statistics
cd("C:/Users/.../PCIJAPTG-A2XvsJSE"); clearconsole()
NPN = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_NPN.csv"); SBK = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_SBK.csv"); SLM = CSV.read("Test Data/A2X/Clean/A2X_Cleaned_SLM.csv")
#---------------------------------------------------------------------------


### 2. Classification rules
function correctmatch(trueclass, inferredclass) # Compare the number of correct classifications against the ground truth. Computes the number of correct classifications.
    n = length(trueclass)
    noclass = 0
    k = 0
    for i in 1:n
        if trueclass[i] == inferredclass[i]
            k += 1
        end
        # Specifically for quote and trade rule. I make the choice s.t.
        # if the algo cannot classify, we do not compare if it can infer
        # therefore reduce the total count of comparisons
        if inferredclass[i] == ""
            noclass += 1
        end
    end
    return k, n, noclass
end
function QuoteRule(data::DataFrame) # obtain the quote rule classification, perform comparisons to obtain percentage of correct classifications
    # Split the data into individual days
    dates = Date.(data[:,2])
    dates_unique = unique(dates)
    # Initialise total counter, correct counter and uncomparable counter
    total = 0
    correct = 0
    uncomparable = 0
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,3])
        # Extract true classifications
        trueclassification = tempdata[tradeinds, 10]
        # Initialise inferred classification
        inferredclassification = String[]
        # Loop through each trade and classify according to quote rule
        for j in 1:length(tradeinds)
            # Get trade value
            TradeValue = tempdata[tradeinds[j], 8]
            tempdata[1:tradeinds[j], 12]
            MidQuoteBeforeTrade = tempdata[findlast(!isnan, tempdata[1:tradeinds[j], 12]), 12]
            # Perform the logic checks. If transaction occurs at a prices higher (lower) than the mid price the trades are classified as buyer (seller) initiated
            if TradeValue > MidQuoteBeforeTrade
                # Transaction is higher than mid price => BuyerInitiated
                push!(inferredclassification, "BuyerInitiated")
            elseif TradeValue < MidQuoteBeforeTrade
                # Transaction is lower than mid price => SellerInitiated
                push!(inferredclassification, "SellerInitiated")
            elseif TradeValue == MidQuoteBeforeTrade
                # Transaction is same as mid quote => no classification
                push!(inferredclassification, "")
            end
        end
        # Compare the classification
        compare = correctmatch(trueclassification, inferredclassification)
        # add the number of comparisons
        total += compare[2]
        # add the number of uncomparable trades
        uncomparable += compare[3]
        # add the number of correct comparisons
        correct += compare[1]
    end
    percentage_correct = correct / (total + uncomparable)
    return percentage_correct, correct, total, uncomparable
end
function TickRule(data::DataFrame) # Obtain the tick rule classification, perform comparisons to obtain percentage of correct classifications
    # Split the data into individual days
    dates = Date.(data[:,2])
    dates_unique = unique(dates)
    # Initialise total counter, correct counter and uncomparable counter
    total = 0
    correct = 0
    uncomparable = 0
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,3])
        # Extract trade values
        TradeValues = tempdata[tradeinds, 8]
        # Extract true classifications
        trueclassification = tempdata[tradeinds, 10]
        # Initialise inferred classification
        inferredclassification = fill("", length(trueclassification), 1)
        # We do not classify the first trade of the day. Don't use the convension of classifying first trade of each day as uptick
        if length(tradeinds) > 1
            # Loop through each trade and classify according to quote rule
            for j in 2:length(tradeinds)
                if TradeValues[j] > TradeValues[j-1]
                    # If trade is higher than previous trade => BuyerInitiated
                    inferredclassification[j] = "BuyerInitiated"
                elseif TradeValues[j] < TradeValues[j-1]
                    # If trade is lower than previous trade => SellerInitiated
                    inferredclassification[j] = "SellerInitiated"

                elseif TradeValues[j] == TradeValues[j-1]
                    # If trades are the same, find last trade that was different
                    indTradeLast = findlast(x -> x != TradeValues[j], TradeValues[1:j])
                    if isnothing(indTradeLast)
                        # No classification, all trades before was the same
                        inferredclassification[j] = ""
                    else
                        # Compare against last trade that was different
                        TradeLast = TradeValues[indTradeLast]
                        if TradeValues[j] > TradeLast
                            # If trade is higher than previous trade => BuyerInitiated
                            inferredclassification[j] = "BuyerInitiated"
                        elseif TradeValues[j] < TradeLast
                            # If trade is lower than previous trade => SellerInitiated
                            inferredclassification[j] = "SellerInitiated"
                        end
                    end
                end
            end
        end
        # Compare the classification
        compare = correctmatch(trueclassification, inferredclassification)
        # Add the number of comparisons
        total += compare[2]
        # Add the number of uncomparable trades
        uncomparable += compare[3]
        # Add the number of correct comparisons
        correct += compare[1]
    end
    percentage_correct = correct / (total + uncomparable)
    return percentage_correct, correct, total, uncomparable
end
function LeeReady(data::DataFrame) # Obtain the Lee-Ready classification, perform comparisons to obtain percentage of correct classifications
    # Split the data into individual days
    dates = Date.(data[:,2])
    dates_unique = unique(dates)
    # Initialise total counter, correct counter and uncomparable counter
    total = 0
    correct = 0
    uncomparable = 0
    @showprogress "Computing..." for i in 1:length(dates_unique)
        # Create extract data from each day
        tempday = dates_unique[i]
        tempdata = data[findall(x -> x == tempday, dates), :]
        # Find the index where trades occur
        tradeinds = findall(x -> x == "TRADE", tempdata[:,3])
        # Extract trade values
        TradeValues = tempdata[tradeinds, 8]
        # Extract true classifications
        trueclassification = tempdata[tradeinds, 10]
        # Initialise inferred classification
        inferredclassification = String[]
        # Loop through each trade and classify according to quote rule
        for j in 1:length(tradeinds)
            # Get trade value
            TradeValue = tempdata[tradeinds[j], 8]
            tempdata[1:tradeinds[j], 12]
            MidQuoteBeforeTrade = tempdata[findlast(!isnan, tempdata[1:tradeinds[j], 12]), 12]
            # Perform the logic checks: begin with quote rule
            if TradeValue > MidQuoteBeforeTrade
                # Transaction is higher than mid price => BuyerInitiated
                push!(inferredclassification, "BuyerInitiated")
            elseif TradeValue < MidQuoteBeforeTrade
                # Transaction is lower than mid price => SellerInitiated
                push!(inferredclassification, "SellerInitiated")
            elseif TradeValue == MidQuoteBeforeTrade
                # Quote rule failed, go to tick rule
                if j > 1
                    if TradeValues[j] > TradeValues[j-1]
                        # If trade is higher than previous trade => BuyerInitiated
                        push!(inferredclassification, "BuyerInitiated")
                    elseif TradeValues[j] < TradeValues[j-1]
                        # If trade is lower than previous trade => SellerInitiated
                        push!(inferredclassification, "SellerInitiated")

                    elseif TradeValues[j] == TradeValues[j-1]
                        # If trades are the same, find last trade that was different
                        indTradeLast = findlast(x -> x != TradeValues[j], TradeValues[1:j])
                        if isnothing(indTradeLast)
                            # No classification, all trades before was the same
                            inferredclassification[j] = ""
                        else
                            # Compare against last trade that was different
                            TradeLast = TradeValues[indTradeLast]
                            if TradeValues[j] > TradeLast
                                # If trade is higher than previous trade => BuyerInitiated
                                inferredclassification[j] = "BuyerInitiated"
                            elseif TradeValues[j] < TradeLast
                                # If trade is lower than previous trade => SellerInitiated
                                inferredclassification[j] = "SellerInitiated"
                            end
                        end
                    end
                else
                    # If first trade can't be classified
                    push!(inferredclassification, "")
                end
            end
        end
        # Compare the classification
        compare = correctmatch(trueclassification, inferredclassification)
        # Add the number of comparisons
        total += compare[2]
        # Add the number of uncomparable trades
        uncomparable += compare[3]
        # Add the number of correct comparisons
        correct += compare[1]
    end
    percentage_correct = correct / (total + uncomparable)
    return percentage_correct, correct, total, uncomparable
end
#---------------------------------------------------------------------------


### 3. Implement classification rules
QuoteRule(NPN); TickRule(NPN); LeeReady(NPN) # Naspers
QuoteRule(SBK); TickRule(SBK); LeeReady(SBK) # Standard Bank
QuoteRule(SLM); TickRule(SLM); LeeReady(SLM) # Sanlam
#---------------------------------------------------------------------------
