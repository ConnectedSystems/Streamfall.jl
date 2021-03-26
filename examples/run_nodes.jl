using YAML, DataFrames, CSV, Plots
using MetaGraphs, Streamfall


# Load and generate stream network
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
g, mg = create_network("Example Network", network)

# Load climate data
climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))
climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."
run_catchment!(mg, g, climate)

inlets, outlets = find_inlets_and_outlets(g)

out = outlets[1]
out_node = get_prop(mg, out, :node)

plot(out_node.outflow)