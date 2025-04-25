using Test
using OrderedCollections

using CSV, YAML
using Dates
using DataFrames
using Streamfall


DATA_PATH = joinpath(dirname(dirname(pathof(Streamfall))), "test/data/hymod")


@testset "Network loading" begin
    # Ensure specified parameter values are being assigned on node creation
    # Load and generate stream network
    sn = load_network("Example Network", joinpath(TEST_DIR, "data/campaspe/campaspe_network.yml"))
    @test sn isa Streamfall.StreamfallNetwork

    # Check known area properties of a few nodes
    msg = """
    Unexpected area value - node may be read in out of order or the network specification
    is not in the expected format.
    """

    target_node = get_prop(sn, 1, :node)
    @test target_node.area == 268.77 || msg

    target_node = get_prop(sn, 2, :node)
    @test target_node.area == 1985.73 || msg

    target_node = get_prop(sn, 8, :node)
    @test target_node.area == 162.84 || msg
end

# Load and generate stream network
sn = load_network("HyMod Network", joinpath(DATA_PATH, "hymod_network.yml"))

@testset "Inlet and outlets of a single node" begin
    inlets, outlets = find_inlets_and_outlets(sn)

    @test length(inlets) == 0 || "No inlets should exist for this specification"
    @test length(outlets) == 1 || "Single node networks should have the node itself as an outlet"
end


@testset "Saving a network spec" begin
    spec = Streamfall.extract_network_spec(sn)
    @test spec isa OrderedDict
    @test haskey(spec, "leaf_river")

    spec["02473000"] = spec["leaf_river"]

    tmp_sn = create_network("Test", spec)

    tmp_fn, io = mktemp()
    Streamfall.save_network(tmp_sn, tmp_fn)
    @test isfile(tmp_fn) || "Network failed to save!"

    # Ensure reloaded network is usable
    new_sn = load_network("Test", tmp_fn)
    @test new_sn isa Streamfall.StreamfallNetwork || "Saved network failed to load"

    # There is currently a bug in YAML.jl which parses string integers as Ints
    # which causes issues when the value starts with a 0
    # (incorrectly parses as an octal integer)
    # If this test does not throw an error, then it means the issue has been resolved
    # and can be removed after checking.
    # see: https://github.com/JuliaData/YAML.jl/pull/45
    @test_throws ErrorException new_sn["02473000"]
end


@testset "Reloading a network spec (HyMod)" begin
    spec::Dict{String,Union{Dict{String,Any},Any}} = Streamfall.extract_network_spec(sn)
    @test spec isa Dict
    @test haskey(spec, "leaf_river") || "Known node does not exist"

    tmp_sn = create_network("Test", spec)

    tmp_fn, io = mktemp()
    Streamfall.save_network(tmp_sn, tmp_fn)
    @test isfile(tmp_fn) || "Network failed to save to file"

    # Ensure reloaded network is usable
    tmp_network = YAML.load_file(tmp_fn)
    new_sn = create_network("Test", tmp_network)

    @test new_sn isa Streamfall.StreamfallNetwork || "Saved network failed to load"

    # Load climate data
    date_format = "YYYY-mm-dd"
    obs_data = CSV.read(
        joinpath(DATA_PATH, "leaf_river_data.csv"),
        DataFrame,
        comment="#",
        dateformat=date_format
    )

    # Column name must match node name
    rename!(obs_data, ["leaf_river_outflow" => "leaf_river"])
    climate_data = obs_data[:, ["Date", "leaf_river_P", "leaf_river_ET"]]
    climate = Climate(climate_data, "_P", "_ET")

    # Ensure reloaded spec still calibrates
    metric = (obs, sim) -> 1.0 - Streamfall.NNSE(obs, sim)
    calibrate!(new_sn, 1, climate, obs_data, metric; MaxTime=5, calibrate_all=false)
end


@testset "Recursing IHACRESNode upstream" begin
    begin
        data_path = joinpath(dirname(dirname(pathof(Streamfall))), "test/data/campaspe/")

        # Load and generate stream network
        sn = load_network("Example Network", joinpath(data_path, "campaspe_network.yml"))

        climate = Climate("../test/data/campaspe/climate/climate.csv", "_rain", "_evap")

        # Historic flows and dam level data
        calib_data = CSV.read(
            joinpath(data_path, "gauges", "outflow_and_level.csv"),
            DataFrame;
            comment="#"
        )

        # Historic extractions from the dam
        extraction_data = CSV.read(
            joinpath(data_path, "gauges", "dam_extraction.csv"),
            DataFrame;
            comment="#"
        )

        @info "Running example stream..."

        dam_id, dam_node = sn["406000"]
        Streamfall.run_node!(sn, dam_id, climate; extraction=extraction_data)

        # Extract data for comparison with 1-year burn-in period
        dam_obs = calib_data[:, "406000"][366:end]
        dam_sim = dam_node.level[366:end]

        nnse_score = Streamfall.NNSE(dam_obs, dam_sim)
        nse_score = Streamfall.NSE(dam_obs, dam_sim)
        rmse_score = Streamfall.RMSE(dam_obs, dam_sim)

        @info "Obj Func Scores:" rmse_score nnse_score nse_score

        nse = round(nse_score, digits=4)
        rmse = round(rmse_score, digits=4)

        # Ensure example does not error out
        reset!(sn)
        run_basin!(sn, climate; extraction=extraction_data)
        @test Streamfall.RMSE(dam_obs, sn[dam_id].level[366:end]) < 2.5


        reset!(sn)
        Streamfall.run_node!(sn, dam_id, climate; extraction=extraction_data)
        @test Streamfall.RMSE(dam_obs, sn[dam_id].level[366:end]) < 2.5

        # Ensure example results haven't changed much...
        # @test nse_score >= 0.95 && nse_score < 1.0 && rmse < 2.0
    end
end