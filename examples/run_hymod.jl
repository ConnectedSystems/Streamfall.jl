using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall


HERE = @__DIR__
DATA_PATH = joinpath(HERE, "../test/data/hymod/")

# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
sn = create_network("HyMod Network", network)

# Load climate data
date_format = "YYYY-mm-dd"
obs_data = CSV.File(joinpath(DATA_PATH, "leaf_river_data.csv"),
                    comment="#",
                    dateformat=date_format) |> DataFrame

hist_streamflow = obs_data[:, ["Date", "leaf_river_outflow"]]
climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
climate = Climate(climate_data, "_P", "_ET")

run_basin!(sn, climate)

nid, node = sn["leaf_river"]
@info mean(node.outflow)