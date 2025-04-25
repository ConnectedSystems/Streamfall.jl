using Test
using Dates, DataFrames, CSV, YAML
using Streamfall


DATA_PATH = joinpath(dirname(dirname(pathof(Streamfall))), "test/data/hymod")

@testset "Single node calibration" begin
    # Load and generate stream network
    network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
    sn = create_network("HyMod Network", network)

    # Load climate data
    date_format = "YYYY-mm-dd"
    obs_data = CSV.File(joinpath(DATA_PATH, "leaf_river_data.csv"),
        comment="#",
        dateformat=date_format) |> DataFrame

    # Column names must be the gauge name
    rename!(obs_data, ["leaf_river_outflow" => "leaf_river"])

    climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
    climate = Climate(climate_data, "_P", "_ET")

    metric = (obs, sim) -> 1.0 - Streamfall.NNSE(obs, sim)

    # `isa Any` check is simply to check that no errors are raised
    # https://github.com/JuliaLang/julia/issues/18780
    # https://discourse.julialang.org/t/test-macro-that-checks-no-exception-was-thrown/41394/2
    @test calibrate!(sn, climate, obs_data, metric; MaxTime=5) isa Any
end
