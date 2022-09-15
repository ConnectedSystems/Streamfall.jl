using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall


here = @__DIR__
data_path = joinpath(here, "../test/data/campaspe/")

# Load and generate stream network
network = YAML.load_file(joinpath(data_path, "campaspe_network.yml"))
sn = create_network("Example Network", network)

# Load climate data
date_format = "YYYY-mm-dd"
climate_data = CSV.File(joinpath(data_path, "climate/climate_historic.csv"),
                          comment="#",
                          dateformat=date_format) |> DataFrame

dam_level_fn = joinpath(data_path, "dam/historic_levels_for_fit.csv")
dam_releases_fn = joinpath(data_path, "dam/historic_releases.csv")
hist_dam_levels = CSV.File(dam_level_fn, dateformat=date_format) |> DataFrame
hist_dam_releases = CSV.File(dam_releases_fn, dateformat=date_format) |> DataFrame

# Combine recorded extractions and releases together
# hist_dam_releases[:, "406000_releases_[ML]"] .+ hist_dam_releases[:, "406000_extractions_[ML]"]
hist_dam_releases[:, "406000_extractions_[ML]"] = hist_dam_releases[:, "406000_releases_[ML]"]

# Subset to same range
climate_data, hist_dam_levels, hist_dam_releases = Streamfall.align_time_frame(climate_data, 
                                                                               hist_dam_levels, 
                                                                               hist_dam_releases)

climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."

reset!(sn)

dam_id, dam_node = sn["406000"]
Streamfall.run_node!(sn, dam_id, climate; extraction=hist_dam_releases)

h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

nnse_score = Streamfall.NNSE(h_data, n_data)
nse_score = Streamfall.NSE(h_data, n_data)
rmse_score = Streamfall.RMSE(h_data, n_data)

@info "Obj Func Scores:" rmse_score nnse_score nse_score

nse = round(nse_score, digits=4)
rmse = round(rmse_score, digits=4)


import Dates: month, monthday, yearmonth

Streamfall.temporal_cross_section(climate_data.Date, h_data, n_data; period=monthday)
savefig("temporal_xsection_monthday_ME.png")

Streamfall.temporal_cross_section(climate_data.Date, h_data, n_data; period=yearmonth)
savefig("temporal_xsection_yearmonth_ME.png")

# Displaying results and saving figure
# plot(h_data,
#      legend=:bottomleft,
#      title="Calibrated IHACRES\n(RMSE: $(rmse); NSE: $(nse))",
#      label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

# display(plot!(n_data, label="IHACRES"))

# savefig("calibrated_example.png")
