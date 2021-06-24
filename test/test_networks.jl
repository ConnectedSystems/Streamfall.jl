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
    spec::Dict{String, Union{Dict{String, Any}, Any}} = Streamfall.extract_network_spec(sn)
    @test spec isa Dict
    @test haskey(spec, "leaf_river")

    spec["02473000"] = spec["leaf_river"]

    tmp_sn = create_network("Test", spec)

    tmp_fn, io = mktemp()
    @test Streamfall.save_network_spec(tmp_sn, tmp_fn) isa Any

    tmp_network = YAML.load_file(tmp_fn)
    tmp_sn = create_network("Test", tmp_network)

    @test tmp_sn isa Streamfall.StreamfallNetwork

    # There is currently a bug in YAML.jl which parses string integers as Ints
    # which causes issues when the value starts with a 0 
    # (incorrectly parses as an octal integer)
    # If this test does not throw an error, then it means the issue has been resolved
    # and can be removed after checking.
    # see:
    # https://github.com/JuliaData/YAML.jl/pull/45
    @test_throws ErrorException sn["02473000"]
end

