# Calibration setup

The calibration examples all rely on the functions shown here.

List of metrics provided by Streamfall can be found in [Included metrics](@ref)


## Importing shared/common packages

```julia
# Ensure dependent data and packages are available
using Statistics, DataFrames, CSV
using Distributed, BlackBoxOptim

using ModelParameters
using LightGraphs, MetaGraphs
using YAML, Plots
using Streamfall
```


## Load network specification

Note that the `DATA_PATH` is pointing to the `test/data/campaspe/` directory.

```julia
# Load and generate stream network
network = YAML.load_file(joinpath(DATA_PATH, "campaspe_network.yml"))
sn = create_network("Example Network", network)
```

## Loading historic data

```julia
# Load climate data
date_format = "YYYY-mm-dd"
climate_data = DataFrame!(CSV.File(joinpath(data_path, "climate/climate_historic.csv"),
                          comment="#",
                          dateformat=date_format))

dam_level_fn = joinpath(data_path, "dam/historic_levels_for_fit.csv")
dam_releases_fn = joinpath(data_path, "dam/historic_releases.csv")
hist_dam_levels = DataFrame!(CSV.File(dam_level_fn, dateformat=date_format))
hist_dam_releases = DataFrame!(CSV.File(dam_releases_fn, dateformat=date_format))

# Subset to same range
climate_data, hist_dam_levels, hist_dam_releases = Streamfall.align_time_frame(climate_data, 
                                                                               hist_dam_levels, 
                                                                               hist_dam_releases)

# Create historic data alias
hist_data = Dict(
    "406000" => hist_dam_levels[:, "Dam Level [mAHD]"]
)

# Create climate object
climate = Climate(climate_data, "_rain", "_evap")
```


## Example objective functions

```julia
"""Calibrate current node."""
function obj_func(params, climate, sn, v_id, calib_data::Dict)

    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    # Running next node will run this node
    Streamfall.run_node!(sn, v_id, climate; extraction=hist_dam_releases)

    n_data = this_node.outflow
    h_data = calib_data[this_node.node_id]

    # Calculate score (NNSE; 0 to 1)
    NNSE = Streamfall.NNSE(h_data, n_data)

    # Switch fitness direction as we want to minimize
    score = 1.0 - NNSE

    # reset to clear stored values
    reset!(sn)

    return score
end


"""Example objective function when performance of current node is dependent 
on the next node.
"""
function obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)

    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    # Run next node which will run this node
    next_node = get_node(sn, next_vid)
    releases = calib_data["$(next_node.node_id)_releases"]
    Streamfall.run_node!(sn, next_vid, climate; extraction=releases)

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


"""Alternative objective function for example. 

This uses a naive split meta-objective function using the Normalized KGE' method.

See `metrics` page for details.
"""
function alt_obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)
    this_node = get_node(sn, v_id)
    update_params!(this_node, params...)

    # Run next node (which will also run this node)
    Streamfall.run_node!(sn, next_vid, climate; extraction=hist_dam_releases)

    next_node = get_node(sn, next_vid)
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
```