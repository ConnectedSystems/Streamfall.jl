using Test
using OrderedCollections

using CSV, YAML
using Dates
using DataFrames
using Streamfall


@testset "Running nodes with interactions" begin
    # Load a network from a file, providing a name for the network and the file path.
    # Creates a graph representation of the stream with associated metadata.
    sn = load_network("Example Network", "../test/data/campaspe/campaspe_network.yml")

    # Load climate data, in this case from a CSV file with data for all nodes.
    climate_data = CSV.read(
        "../test/data/campaspe/climate/climate_historic.csv",
        DataFrame;
        comment="#",
        dateformat="YYYY-mm-dd"
    )

    # Indicate which columns are precipitation and evaporation data based on partial identifiers
    climate = Climate(climate_data, "_rain", "_evap")

    run_catchment!(sn, climate)
    node_id, node = sn["406219"]

    # Run up to a point in the stream for all time steps.
    # All nodes upstream will be run as well (but not those downstream)
    run_node!(sn, node_id, climate)

    baseline_outflow = node.outflow

    # Reset a node (clears stored states)
    reset!(node)

    # Run a specific node, and only a specific node, for all time steps
    inflow = 10.0      # inflows for each time step
    extractions = 0.0 # extractions from stream for each time step
    gw_flux = 0.0     # forced groundwater interactions for each time step
    run_node!(node, climate; inflow=inflow, extraction=extractions, exchange=gw_flux)
    perturbed1 = node.outflow

    reset!(node)
    inflow = 10.0      # inflows for each time step
    extractions = 5.0 # extractions from stream for each time step
    gw_flux = 0.0     # forced groundwater interactions for each time step
    run_node!(node, climate; inflow=inflow, extraction=extractions, exchange=gw_flux)
    perturbed2 = node.outflow

    reset!(sn)
    inflow = 10.0      # inflows for each time step
    extractions = 5.0 # extractions from stream for each time step
    gw_flux = 2.0     # forced groundwater interactions for each time step
    inlets, outlets = find_inlets_and_outlets(sn)

    # Manual interactions
    timesteps = sim_length(climate)
    prep_state!(sn, timesteps)
    for ts in (1:timesteps)
        for outlet in outlets
            run_node!(sn[outlet], climate, ts; inflow=inflow, extraction=extractions, exchange=gw_flux)
        end
    end

    node_id, node = sn["406219"]
    perturbed3 = node.outflow

    @test all(baseline_outflow .!= perturbed1 .!= perturbed2 .!= perturbed3) || "Perturbations resulted in identical streamflows"

    # Check that additional interactions are accounted for.
    # Cannot check for equivalence due to computational error, so we test that the
    # difference is within acceptable bounds.
    @test all((baseline_outflow .+ 10.0 .- perturbed1) .< 0.0005) || "Additional inflow not accounted for"
    @test all((baseline_outflow .+ 5.0 .- perturbed2) .< 0.0005) || "Extractions not accounted for"
end

