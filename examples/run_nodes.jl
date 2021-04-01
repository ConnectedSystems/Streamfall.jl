using YAML, DataFrames, CSV, Plots
using MetaGraphs, Streamfall


# Load and generate stream network
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
g, mg = create_network("Example Network", network)

# Load climate data
climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_fortran.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))
climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."
run_catchment!(mg, g, climate)

inlets, outlets = find_inlets_and_outlets(g)

in_node = get_prop(mg, inlets[1], :node)

plot(in_node.outflow)

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
