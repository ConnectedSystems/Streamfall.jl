using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall


here = @__DIR__
data_path = joinpath(here, "../test/data/campaspe/")

# Load and generate stream network
sn = load_network("Example Network", joinpath(data_path, "campaspe_network.yml"))

# The Campaspe catchment is represented as a network of eight nodes, including one dam.
# All nodes use the IHACRES_CMD rainfall-runoff model.
plot_network(sn)

# Load climate data - in this case from a CSV file with data for all nodes.
climate_data = CSV.read(
    joinpath(data_path, "climate", "climate.csv"),
    DataFrame;
    comment="#"
)

# Indicate which columns are precipitation and evaporation data based on partial identifiers
climate = Climate(climate_data, "_rain", "_evap")

# Historic flows and dam level data
calib_data = CSV.read(
    joinpath(data_path, "gauges", "outflow_and_level.csv"),
    DataFrame;
    comment="#"
)

# Historic extractions from the dam
extraction_data = CSV.read(
    joinpath(data_path, "gauges", "dam_extraction.csv"),
    DataFrame;
    comment="#"
)

@info "Running example stream..."

dam_id, dam_node = sn["406000"]
Streamfall.run_node!(sn, dam_id, climate; extraction=extraction_data)

# Extract data for comparison with 1-year burn-in period
dam_obs = calib_data[:, "406000"][366:end]
dam_sim = dam_node.level[366:end]

nnse_score = Streamfall.NNSE(dam_obs, dam_sim)
nse_score = Streamfall.NSE(dam_obs, dam_sim)
rmse_score = Streamfall.RMSE(dam_obs, dam_sim)

@info "Obj Func Scores:" rmse_score nnse_score nse_score

nse = round(nse_score, digits=4)
rmse = round(rmse_score, digits=4)

import Dates: month, monthday, yearmonth

climate_burnin = climate_data[366:end, :Date]
Streamfall.temporal_cross_section(climate_burnin, dam_obs, dam_sim; period=monthday)
# savefig("temporal_xsection_monthday_ME.png")

Streamfall.temporal_cross_section(climate_burnin, dam_obs, dam_sim; period=yearmonth)
# savefig("temporal_xsection_yearmonth_ME.png")

# Displaying results and saving figure
# plot(dam_obs,
#      legend=:bottomleft,
#      title="Calibrated IHACRES\n(RMSE: $(rmse); NSE: $(nse))",
#      label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

# display(plot!(dam_sim, label="IHACRES"))

# savefig("calibrated_example.png")
