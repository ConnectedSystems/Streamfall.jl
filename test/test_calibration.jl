using Test
using Dates, DataFrames, CSV, YAML
using Streamfall


DATA_PATH = joinpath(@__DIR__, "data/hymod") 

# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
sn = create_network("HyMod Network", network)

# Load climate data
date_format = "YYYY-mm-dd"
obs_data = DataFrame!(CSV.File(joinpath(DATA_PATH, "leaf_river_data.csv"),
                          comment="#",
                          dateformat=date_format))

hist_streamflow = obs_data[:, ["Date", "leaf_river_outflow"]]
climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
climate = Climate(climate_data, "_P", "_ET")


@testset "Single node calibration" begin
    @test calibrate!(sn, climate, hist_streamflow; MaxTime=30) isa Any
end
