using YAML
using DataFrames, CSV

using Streamfall


HERE = @__DIR__

include(HERE*"/network_creation.jl")

climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))
climate = Climate(climate_data, "_rain", "_evap")

@info "Running example stream..."
timesteps = sim_length(climate)
for ts in (1:timesteps)
    for outlet in outlets
        run_node!(mg, g, outlet, climate, ts)
    end
end

using MetaGraphs, Plots

out = outlets[1]
out_node = get_prop(mg, out, :node)



plot(out_node.outflow)
