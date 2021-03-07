using YAML
using LightGraphs, MetaGraphs
using DataFrames, CSV

import ModelParameters: update, update!, Model
using Streamfall

using Infiltrator


network = YAML.load_file("../tests/data/AWRA_R_Network/campaspe_network.yml")
g, mg = create_network("Example Network", network)

inlets, outlets = find_inlets_and_outlets(g)

@info "Network has the following inlets and outlets:" inlets outlets



# climate_data = DataFrame!(CSV.File("../tests/data/climate/climate_historic.csv", 
#                           comment="#",
#                           dateformat="YYYY-mm-dd"))
# climate = Climate(climate_data, "_rain", "_evap")

# @info "Running example stream..."
# timesteps = sim_length(climate)
# for ts in (1:timesteps)
#     for outlet in outlets
#         run_node!(mg, g, outlet, climate, ts)
#     end
# end