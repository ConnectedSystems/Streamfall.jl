using Test
using Dates, DataFrames, CSV, YAML
using Streamfall


DATA_PATH = joinpath(@__DIR__, "data/hymod") 

# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "hymod_network.yml"))
sn = create_network("HyMod Network", network)


@testset "Inlet and outlets" begin
    inlets, outlets = find_inlets_and_outlets(sn)

    @test length(inlets) == 0
    @test length(outlets) == 1
end