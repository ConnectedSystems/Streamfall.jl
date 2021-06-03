using YAML, DataFrames, CSV, Plots
using Statistics
using MetaGraphs, Streamfall


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

dam_id, dam_node = get_gauge(mg, "406000")
timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(sn, dam_id, climate, ts; water_order=hist_dam_releases)
end

h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

@info "NNSE:" Streamfall.NNSE(h_data, n_data)
@info "NSE:" Streamfall.NSE(h_data, n_data)
@info "RMSE:" Streamfall.RMSE(h_data, n_data)

nse = round(Streamfall.NSE(h_data, n_data), digits=4)
rmse = round(Streamfall.RMSE(h_data, n_data), digits=4)

plot(h_data,
     legend=:bottomleft,
     title="IHACRES Calibration\n(NSE: $(nse); RMSE: $(rmse))",
     label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

plot!(n_data, label="IHACRES")

savefig("calibration_ts_comparison.png")

# 1:1 Plot
scatter(h_data, n_data, legend=false, 
        markerstrokewidth=0, markerstrokealpha=0, alpha=0.2)
plot!(h_data, h_data, color=:red, markersize=.1, markerstrokewidth=0,
      xlabel="Historic [mAHD]", ylabel="IHACRES [mAHD]", title="Historic vs Modelled")

savefig("calibration_1to1.png")
