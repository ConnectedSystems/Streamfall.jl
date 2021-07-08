using Test
using Dates, DataFrames, CSV, YAML
using Streamfall


DATA_PATH = joinpath(@__DIR__, "data/hymod") 

# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
sn = create_network("HyMod Network", network)


@testset "Single node network" begin
    inlets, outlets = find_inlets_and_outlets(sn)

    @test length(inlets) == 0
    @test length(outlets) == 1
end


@testset "Saving a network spec" begin
    spec = Streamfall.extract_network_spec(sn)
    @test spec isa Dict
    @test haskey(spec, "leaf_river")

    spec["02473000"] = spec["leaf_river"]

    tmp_sn = create_network("Test", spec)

    tmp_fn, io = mktemp()
    @test Streamfall.save_network_spec(tmp_sn, tmp_fn) isa Any

    # Ensure reloaded network is usable
    tmp_network = YAML.load_file(tmp_fn)
    new_sn = create_network("Test", tmp_network)

    @test new_sn isa Streamfall.StreamfallNetwork

    # There is currently a bug in YAML.jl which parses string integers as Ints
    # which causes issues when the value starts with a 0 
    # (incorrectly parses as an octal integer)
    # If this test does not throw an error, then it means the issue has been resolved
    # and can be removed after checking.
    # see: https://github.com/JuliaData/YAML.jl/pull/45
    @test_throws ErrorException new_sn["02473000"]
end


@testset "Reloading a network spec" begin
    spec::Dict{String, Union{Dict{String, Any}, Any}} = Streamfall.extract_network_spec(sn)
    @test spec isa Dict
    @test haskey(spec, "leaf_river")

    tmp_sn = create_network("Test", spec)

    tmp_fn, io = mktemp()
    @test Streamfall.save_network_spec(tmp_sn, tmp_fn) isa Any

    # Ensure reloaded network is usable
    tmp_network = YAML.load_file(tmp_fn)
    new_sn = create_network("Test", tmp_network)

    @test new_sn isa Streamfall.StreamfallNetwork

    # Load climate data
    date_format = "YYYY-mm-dd"
    obs_data = DataFrame!(CSV.File(joinpath(DATA_PATH, "leaf_river_data.csv"),
                            comment="#",
                            dateformat=date_format))

    hist_streamflow = obs_data[:, ["Date", "leaf_river_outflow"]]
    climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
    climate = Climate(climate_data, "_P", "_ET")

    # Ensure reloaded spec still calibrates
    metric = (obs, sim) -> 1.0 - Streamfall.NNSE(obs, sim)
    @test calibrate!(new_sn, climate, hist_streamflow; metric=metric, MaxTime=10) isa Any
end


@testset "Recursing IHACRESNode upstream" begin

    begin
        include("../examples/run_nodes.jl")
        # Ensure other methods of running a node are identical

        reset!(sn)
        run_basin!(sn, climate; extraction=hist_dam_releases)
        @test Streamfall.RMSE(n_data, sn[dam_id].level) == 0.0


        reset!(sn)
        Streamfall.run_node!(sn, dam_id, climate; extraction=hist_dam_releases)
        @test Streamfall.RMSE(n_data, sn[dam_id].level) == 0.0

        @test rmse >= 0.95
    end
end