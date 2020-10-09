## Author: Patrick Chang
# Script file to read in the A2X trade data, extract
# the usable messeages, write data into clean dataframes,
# and write the data into flat files

## Phase 1 of cleaning
#---------------------------------------------------------------------------
## Preamble
using CSV, DataTables, CodecBzip2, DataFrames, JLD, Dates, ProgressMeter

cd("/Users/patrickchang1/HFT2020")

# Obtain Dictionary of securityId to the name of the ticker
# CSV.write("SecurityIDtoTickerName.csv", secIDtoTickerName)
SecurityIDtoTickerName = CSV.read("SecurityIDtoTickerName.csv")
secIDtoTickerName = Dict()
for i in 1:size(SecurityIDtoTickerName)[1]
    push!(secIDtoTickerName, SecurityIDtoTickerName[i,1] => SecurityIDtoTickerName[i,2])
end

#---------------------------------------------------------------------------
## Functions to clean the data
# Function to extract value of a specific field from a message string
function getFromString(field::String, message::SubString{String}, delimiter::String)
    first_ind = findfirst(field*":", message)[end]+1
    last_ind = findfirst(delimiter, message[first_ind:end])[1]-2
    return message[first_ind:first_ind+last_ind]
end
# Function to convert unix time with nanoseconds into DateTime with milliseconds
function Unix2DateTime(time::Float64)
    seconds = time ÷ 10^9
    milliseconds = Millisecond(time ÷ 10^6 - seconds * 10^3)
    nanoseconds = Nanosecond(time % 10^6)
    return unix2datetime(seconds + 2*3600) + milliseconds + nanoseconds
end
# Function to create dataframe information from a message string.
# The structure of fields to enter depend on the type of message string.
function makeDataStruct(message::SubString{String}, type)
    if type == 2
        securityId = getFromString("securityId", message, ",")
        side = getFromString("side", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("limitPrice", message, ","))
        orderId = getFromString("orderId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OA", price, quantity, side, "", "")
    elseif type == 3
        securityId = getFromString("securityId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        orderId = getFromString("orderId", message, ",")
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OC", NaN, NaN, "", "", "")
    elseif type == 4
        securityId = getFromString("securityId", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        orderId = getFromString("orderId", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "OM", price, quantity, "", "", "")
    elseif type == 5
        securityId = getFromString("securityId", message, ",")
        tradeType = getFromString("tradeType", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        orderId = getFromString("orderId", message, ",")
        tradeRef = getFromString("tradeRef", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, orderId, "AT", price, quantity, "", tradeRef, tradeType)
    else    # type 6 message
        securityId = getFromString("securityId", message, ",")
        quantity = parse(Float64, getFromString("quantity", message, ","))
        price = parse(Float64, getFromString("price", message, ","))
        tradeRef = getFromString("tradeRef", message, ",")
        timestamp = parse(Float64, getFromString("timestamp", message, "}"))
        date = Unix2DateTime(timestamp)
        return (securityId, timestamp, date, "", "TB", price, NaN, "", tradeRef, "")
    end
end
# Function to read in a day of data of all assets, remove useless information
# construct dataframe of useful information for each asset. Each DataFrame is
# stored as a value in a dictionary.
function getTAQ(file::String)
    ## Load in the Data
    compressed = read(file)
    decompress = transcode(Bzip2Decompressor, compressed)
    data = String(decompress)

    data_split = split(data, "\n")

    ## Two sets of messeages to:
    MdHeader = "MdHeader"   # Keep
    type1 = "msgType:1"     # Remove

    # Kepp all MdHeader messages
    indkeep = occursin.(MdHeader, data_split)
    data_split = data_split[indkeep]

    # Remove all Type 1 messages
    indkeep2 = occursin.(type1, data_split)
    indkeep2 = .!indkeep2
    data_split = data_split[indkeep2]

    data_split2 = []
    for i in 1:length(data_split)
        push!(data_split2, split(data_split[i], " ")[2])
    end

    ## Initialise DataFrames inside dictionary
    type8 = "msgType:8"
    indkeepT8 = occursin.(type8, data_split2)
    mktInfo = data_split2[indkeepT8]
    mktInfo_split = split.(mktInfo,"{")

    # Get securityId and Name of ticker
    secID = []
    tickerName = []
    for i in 1:length(mktInfo_split)
        push!(secID, split(mktInfo_split[i][4], ",")[1])
        push!(tickerName, split(mktInfo_split[i][4], ",")[2][7:end-2])
    end

    # # Create dictionary mapping form securityId to tickerName
    # secIDtoTickerName = Dict()
    # for i in 1:length(secID)
    #     push!(secIDtoTickerName, secID[i] => tickerName[i])
    # end

    # Data Structure for each row of DataFrame for all assets
    df = DataFrame(Security = String[], TimeStamp = Float64[], Date = DateTime[],
    OrderRef = String[],Type = String[], Price = Float64[], Quantity = Float64[],
    Side = String[], TradeRef = String[], TradeType = String[])

    # Loop through data set
    for i in 1:length(data_split2)
        # Identify securityId and msgType
        message = data_split2[i]
        msgType = parse(Int, getFromString("msgType", message, ","))

        # if appropriate message
        if (msgType == 2 || msgType == 3 || msgType == 4 || msgType == 5 || msgType == 6)
            # Extract DataFrame format to push into dataframe
            getDataStruct = makeDataStruct(message, msgType)
            push!(df, getDataStruct)
        end
    end

    # Create dictionary of DataFrames for each asset
    df_store = Dict()
    for i in 1:length(secID)
        Id = split(secID[i], ":")[2]
        ind_Id = findall(x -> x == Id, df[:Security])
        if length(ind_Id) == 0
            df_empty = DataFrame(Security = String[], TimeStamp = Float64[], Date = DateTime[],
            OrderRef = String[],Type = String[], Price = Float64[], Quantity = Float64[],
            Side = String[], TradeRef = String[], TradeType = String[])
            push!(df_store, secIDtoTickerName[secID[i]] => df_empty)
        else
            df_temp = df[ind_Id,:]
            push!(df_store, secIDtoTickerName[secID[i]] => df_temp)
        end
    end
    return df_store, secID, tickerName
end
# Function to streamline everything. Reads in all the data and returns a dictionary
# of dataframes for the entire history of the A2X equities.
function combineTAQ(filename)
    Master_df = Dict()
    @showprogress "Computing..." for i in 1:length(filename)
        day_data = getTAQ(filename[i])
        df_day = day_data[1]
        ticker = day_data[3]

        for j in 1:length(ticker)
            if haskey(Master_df, ticker[j])
                temp = [ Master_df[ticker[j]]; df_day[ticker[j]] ]
                Master_df[ticker[j]] = temp
            else
                push!(Master_df, ticker[j] => df_day[ticker[j]])
            end
        end
    end
    return Master_df
end

#---------------------------------------------------------------------------
## Extract useful information from string of data
files = readdir("A2X/Raw/")
filename = "A2X/Raw/".*files
filename = filename     # should be 911 items, if more remove DSstore.
# Roughly takes 1 day and 1hour 30 minute to compute
Full_data = combineTAQ(filename)

#---------------------------------------------------------------------------
## Write the cleaned dataframes into flat files
TickerName = []
for i in 1:size(SecurityIDtoTickerName)[1]
    push!(TickerName, SecurityIDtoTickerName[i,2])
end

@showprogress "Computing..." for i in 1:length(TickerName)
    CSV.write("A2X/TAQ_full/TAQ_FULL_"*TickerName[i]*".csv", Full_data[TickerName[i]])
end
# Save full data in jld format
save("A2X/TAQ_full/Full_data.jld", Full_data)
