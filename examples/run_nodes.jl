using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall


# Load and generate stream network
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
sn = create_network("Example Network", network)

# Load climate data
climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv",
                          comment="#",
                          dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

# Subset to same range
first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1], climate_data.Date[1])
last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end], climate_data.Date[end])

climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]

climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."

reset!(sn)

dam_id, dam_node = get_gauge(sn, "406000")
Streamfall.run_node!(sn, dam_id, climate; water_order=hist_dam_releases)

h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

nnse_score = Streamfall.NNSE(h_data, n_data)
nse_score = Streamfall.NSE(h_data, n_data)
rmse_score = Streamfall.RMSE(h_data, n_data)

@info "NNSE:" nnse_score
@info "NSE:" nse_score
@info "RMSE:" rmse_score

nse = round(nse_score, digits=4)
rmse = round(rmse_score, digits=4)

plot(h_data,
     legend=:bottomleft,
     title="Calibrated IHACRES\n(NSE: $(nse); RMSE: $(rmse))",
     label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

plot!(n_data, label="IHACRES")
