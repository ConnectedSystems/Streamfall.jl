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

    tmp_fn, io = mktemp()
    @test Streamfall.save_network_spec(sn, tmp_fn) isa Any

    tmp_network = YAML.load_file(tmp_fn)
    tmp_sn = create_network("Test", tmp_network)

    @test tmp_sn isa Streamfall.StreamfallNetwork

    @info sn["leaf_river"]
    @info tmp_sn["leaf_river"]
end

