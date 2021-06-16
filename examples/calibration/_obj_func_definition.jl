DATA_PATH = "../../test/data/campaspe/"

# Ensure dependent data and packages are available
using Statistics, DataFrames, CSV
using Distributed, BlackBoxOptim

using ModelParameters
using LightGraphs, MetaGraphs
using YAML
using Streamfall

import Random


Random.seed!(101)
network = YAML.load_file(joinpath(DATA_PATH, "campaspe_network.yml"))
sn = create_network("Example Network", network)

# inlets, outlets = find_inlets_and_outlets(sn)
# @info "Network has the following inlets and outlets:" inlets outlets

climate_data = DataFrame!(CSV.File(joinpath(DATA_PATH, "climate/climate_historic.csv"),
                            comment="#",
                            dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File(joinpath(DATA_PATH, "dam/historic_levels_for_fit.csv"), dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File(joinpath(DATA_PATH, "dam/historic_releases.csv"), dateformat="YYYY-mm-dd"))

# We're ignoring flow data as it is severely limited (only ~5 years)
# hist_406219_flow = DataFrame!(CSV.File(joinpath(DATA_PATH, "gauges/406219_outflow_edited.csv"), dateformat="YYYY-mm-dd"))

# Subset climate and historic flow data to same range
first_date, last_date = Streamfall.find_common_timeframe(hist_dam_levels, hist_dam_releases, climate_data)

climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]
hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]
# hist_406219_flow = hist_406219_flow[first_date .<= hist_406219_flow.Date .<= last_date, :]

hist_data = Dict(
    # "406219" => hist_406219_flow[:, "406219_outflow_[ML]"],
    "406000" => hist_dam_levels[:, "Dam Level [mAHD]"],
    "406000_releases" => hist_dam_releases
)

climate = Climate(climate_data, "_rain", "_evap")


"""Calibrate current node."""
function obj_func(params, climate, sn, v_id, calib_data::Dict)

    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    releases = calib_data["$(this_node.node_id)_releases"]
    Streamfall.run_node!(sn, v_id, climate; water_order=releases)

    n_data = this_node.outflow
    h_data = calib_data[this_node.node_id]

    # Calculate score (NNSE; 0 to 1)
    NNSE = Streamfall.NNSE(h_data, n_data)

    # Switch fitness direction as we want to minimize
    score = 1.0 - NNSE

    # Other metrics are also available.
    # see `metrics.jl`
    # RMSE = Streamfall.RMSE(h_data, node_data)
    # score = RMSE

    # reset to clear stored values
    reset!(sn)

    return score
end


"""Example objective function when performance of current node is
dependent on the next node.
"""
function obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)

    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    # Run next node which will run this node
    next_node = get_node(sn, next_vid)
    releases = calib_data["$(next_node.node_id)_releases"]
    Streamfall.run_node!(sn, next_vid, climate; water_order=releases)

    # Alias data as necessary
    if next_node.node_id == "406000"
        n_data = next_node.level
        h_data = calib_data[next_node.node_id]
    elseif this_node.node_id == "406000"
        n_data = this_node.level
        h_data = calib_data[this_node.node_id]
    else
        n_data = this_node.outflow
        h_data = calib_data[this_node.node_id]
    end

    NNSE = Streamfall.NNSE(h_data, n_data)
    score = 1.0 - NNSE

    reset!(sn)

    return score
end


function alt_obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)
    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    next_node = get_node(sn, next_vid)

    # Run next node (which will also run this node)
    releases = calib_data["$(next_node.node_id)_releases"]
    Streamfall.run_node!(sn, next_vid, climate; water_order=releases)

    # Alias data as necessary
    if next_node.node_id == "406000"
        n_data = next_node.level
        h_data = calib_data[next_node.node_id]
    elseif this_node.node_id == "406000"
        n_data = this_node.level
        h_data = calib_data[this_node.node_id]
    else
        n_data = this_node.outflow
        h_data = calib_data[this_node.node_id]
    end

    split_NmKGE = Streamfall.naive_split_metric(h_data, n_data; n_members=365, metric=Streamfall.NmKGE, comb_method=mean)
    score = 1.0 - split_NmKGE

    reset!(sn)

    return score
end