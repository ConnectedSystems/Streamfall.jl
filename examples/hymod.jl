using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall

here = @__DIR__
data_path = joinpath(here, "../test/data/hymod/")

# Load and generate stream network
network = YAML.load_file(joinpath(data_path, "hymod_network.yml"))
sn = create_network("HyMod Network", network)

# Load climate data
date_format = "YYYY-mm-dd"
obs_data = DataFrame!(CSV.File(joinpath(data_path, "leaf_river_data.csv"),
                          comment="#",
                          dateformat=date_format))

hist_streamflow = obs_data[:, ["Date", "leaf_river_streamflow"]]
climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
climate = Climate(climate_data, "_P", "_ET")

run_basin!(sn, climate)

nid, node = get_gauge(sn, "leaf_river")
@info mean(node.outflow)