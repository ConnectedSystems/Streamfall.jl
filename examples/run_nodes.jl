using YAML, DataFrames, CSV, Plots
using Statistics
using MetaGraphs, Streamfall


# Load and generate stream network
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
mg, g = create_network("Example Network", network)

# Load climate data
climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

# Subset to same range
first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1])
last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end])

climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]

climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."
run_catchment!(mg, g, climate; water_order=hist_dam_releases)

match = collect(filter_vertices(mg, :name, "406000"))
dam_id = match[1]
dam_node = get_prop(mg, dam_id, :node)

h_levels = hist_dam_levels[:, "Dam Level [mAHD]"]
n_levels = dam_node.level

@info "NNSE:" Streamfall.NNSE(h_levels, n_levels)
@info "RMSE:" Streamfall.RMSE(h_levels, n_levels)


plot(h_levels)
display(plot!(n_levels))

# outflow = in_node.outflow
# append!(outflow, NaN)
# erain = in_node.effective_rainfall
# append!(erain, NaN)

# res = Dict(
#     "Outflow"=> outflow,
#     "CMD"=>in_node.storage,
#     "Quick Store"=>in_node.quick_store,
#     "Slow Store"=>in_node.slow_store,
#     "erain"=>erain
# )

# CSV.write("outflow.csv", DataFrame(res))
