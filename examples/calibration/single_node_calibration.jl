using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall


DATA_PATH = joinpath(dirname(dirname(pathof(Streamfall))), "test/data/hymod/")

# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
sn = create_network("HyMod Network", network)

# Load climate data
obs_data = CSV.read(
    joinpath(DATA_PATH, "leaf_river_data.csv"),
    DataFrame,
    comment="#",
    dateformat="YYYY-mm-dd"
)

hist_streamflow = obs_data[:, "leaf_river_outflow"]
rename!(obs_data, ["leaf_river_outflow" => "leaf_river"])

climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
climate = Climate(climate_data, "_P", "_ET")

# This will set node parameters to the optimal values found
metric = (obs, sim) -> 1.0 - Streamfall.NNSE(obs, sim)
calibrate!(sn, climate, obs_data, metric; MaxTime=180.0)

# Save calibrated network spec to file
Streamfall.save_network(sn, "hymod_example_calibrated.yml")

# Run the model as a single system
run_basin!(sn, climate)

# Get node
nid, node = sn["leaf_river"]

@info "Mean flow (ft^3/s)" mean(node.outflow)

# Get performance and plot
obs = hist_streamflow[:, "leaf_river_outflow"]
@info "RMSE" Streamfall.RMSE(obs, node.outflow)
@info "NSE" Streamfall.NSE(obs, node.outflow)

plot(obs)
plot!(node.outflow)
