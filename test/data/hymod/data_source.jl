# Download data for USGS 02473000 - Leaf River at Hattiesburg, MS
# variable=00060,00065 refers to discharge in cubic feet per second and gage height in feet.

using JSON, CSV, DataFrames

# Request discharge and gage height data
json_data = download("https://waterservices.usgs.gov/nwis/dv/?site=02473000&startDT=1981-01-01&endDT=2021-06-22&variable=00060,00065&format=json")

data = open(json_data, "r")

parsed = JSON.Parser.parse(data)

time_series = parsed["value"]["timeSeries"]

# `time_series` is an array, where the first entry is discharge/streamflow and the second is gage height.

# Metadata for streamflow
# time_series[1]["variable"]

# time_series[1]["values"][1]
# actual time series
streamflow_raw = time_series[1]["values"][1]["value"]
level_raw = time_series[2]["values"][1]["value"]

streamflow_ft = vcat(DataFrame.(streamflow_raw)...)
level_ft = vcat(DataFrame.(level_raw)...)

stream_col = Symbol("leaf_river_outflow_[ft^3/s]")
level_col = Symbol("leaf_river_level_[ft]")
rename!(streamflow_ft, :value => stream_col)
rename!(level_ft, :value => level_col)

# Combine the time series
timeseries = innerjoin(streamflow_ft[:, [:dateTime, stream_col]], level_ft[:, [:dateTime, level_col]], on=:dateTime)

# Strip time info from date
timeseries.dateTime = getindex.(string.(timeseries.dateTime), Ref(1:10))
rename!(timeseries, :dateTime => :Date)


CSV.write("leaf_river.csv", timeseries)
