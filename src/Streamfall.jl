module Streamfall

using Statistics
using Graphs, MetaGraphs, Distributed, DataFrames


const MODPATH = @__DIR__

if Sys.iswindows()
    libext = ".dll"
elseif Sys.islinux()
    libext = ".so"
elseif Sys.isapple()
    libext = ".dynlib"
else
    throw(DomainError("Unsupported platform"))
end

# Can't use string, DLL location has to be a const
# (which makes sense but still, many hours wasted!)
# https://github.com/JuliaLang/julia/issues/29602
const IHACRES = joinpath(MODPATH, "../deps", "usr", "lib", "ihacres$(libext)")


include("Network.jl")
include("Nodes/Node.jl")
include("Climate.jl")
include("metrics.jl")
include("calibration.jl")
include("Nodes/IHACRES/IHACRESNode.jl")
include("Nodes/IHACRES/IHACRESExpuhNode.jl")
include("Nodes/GR4J/GR4JNode.jl")
include("Nodes/HyMod/HyModNode.jl")
include("Nodes/SYMHYD/SYMHYDNode.jl")
include("Nodes/DamNode.jl")
include("Nodes/Ensembles/EnsembleNode.jl")


function timestep_value(ts::Int, gauge_id::String, col_partial::String, dataset::DataFrame)::Float64
    target_col = filter(x -> occursin(gauge_id, x)
                                & occursin(col_partial, x),
                                names(dataset)
    )

    amount = 0.0
    if !isempty(target_col)
        amount = if checkbounds(Bool, dataset.Date, ts)
            dataset[ts, target_col][1]
        else
            0.0
        end
    end

    return amount
end
@inline function timestep_value(_::Int, _::String, _::String, dataset::Float64)::Float64
    return dataset
end
@inline function timestep_value(_::Int, _::String, _::String, dataset::Nothing)::Float64
    return 0.0
end
@inline function timestep_value(ts::Int, _::String, _::String, dataset::Vector)::Float64
    return dataset[ts]
end


"""
    find_common_timeframe(timeseries::T...)

Find common time frame between time series.

Requires that all DataFrames have a "Date" column.
"""
function find_common_timeframe(timeseries::T...) where {T<:DataFrame}
    for t in timeseries
        @assert "Date" in names(t)
    end

    min_date = maximum([x.Date[1] for x in timeseries])
    max_date = minimum([x.Date[end] for x in timeseries])
    return (min_date, max_date)
end

"""
    prep_state!(sn::StreamfallNetwork, timesteps::Int64)::Nothing

Prepare a network for a run by pre-allocating result stores.
"""
function prep_state!(sn::StreamfallNetwork, timesteps::Int64)::Nothing
    for node in sn
        prep_state!(node, timesteps)
    end
end


"""
    align_time_frame(timeseries::T...)

Subset an arbitrary number of DataFrames to their shared period of time.

Returns subsetted copy of data in same order as input.

# Example
```julia-repl
julia> climate, streamflow = align_time_frame(climate, streamflow)
```
"""
function align_time_frame(timeseries::T...) where {T<:DataFrame}
    min_date, max_date = find_common_timeframe(timeseries...)
    modded = [t[min_date .<= t.Date .<= max_date, :] for t in timeseries]

    for t in modded
        ts_diff = diff(t.Date)
        is_contiguous = all(ts_diff .== ts_diff[1])

        if !is_contiguous
            msg = """
            Data is non-contiguous or is inconsistent (e.g., a mix of daily and monthly data).
            There cannot be gaps in the time series.
            """
            throw(ArgumentError(msg))
        end
    end

    return modded
end

"""
    run_basin!(sn::StreamfallNetwork, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)

Run scenario for an entire catchment/basin.
"""
function run_basin!(
    sn::StreamfallNetwork, climate::Climate;
    inflow=nothing, extraction=nothing, exchange=nothing
)
    _, outlets = find_inlets_and_outlets(sn)
    @inbounds for outlet_id in outlets
        run_node!(sn, outlet_id, climate;
                  inflow=inflow, extraction=extraction, exchange=exchange)
    end
end

run_catchment! = run_basin!


"""
    run_node!(
        sn::StreamfallNetwork, node_id::Int, climate::Climate, ts::Int64;
        inflow=nothing, extraction=nothing, exchange=nothing
    )::Nothing

Generic run method that runs a model attached to a given node for a given timestep.
Recurses upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id` : Node to run in the network
- `climate` : Climate object holding rainfall and evaporation data (or temperature)
- `ts` : Timestep to run
- `extraction` : Water orders for each time step (defaults to nothing)
- `exchange` : Exchange with groundwater system at each time step (defaults to nothing)
"""
function run_node!(
    sn::StreamfallNetwork, node_id::Int, climate::Climate, ts::Int64;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Nothing
    timesteps = sim_length(climate)

    # Run all upstream nodes
    sim_inflow = zeros(timesteps)
    ins = inneighbors(sn.mg, node_id)
    for i in ins
        # Get inflow from previous node
        run_node!(
            sn, i, climate, ts;
            inflow=inflow, extraction=extraction, exchange=exchange
        )

        # Add outflow from upstream to inflow
        sim_inflow .+= sn[i].outflow
    end

    # Run this node
    run_node!(
        sn[node_id], climate, ts;
        inflow=sim_inflow, extraction=extraction, exchange=exchange
    )

    return nothing
end

function run_node!(
    sn::StreamfallNetwork, node_id::Int, climate::Climate;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Nothing
    timesteps = sim_length(climate)

    # Run all upstream nodes
    sim_inflow = zeros(timesteps)
    ins = inneighbors(sn.mg, node_id)
    for i in ins
        # Get inflow from previous node
        run_node!(
            sn, i, climate;
            inflow=inflow, extraction=extraction, exchange=exchange
        )

        # Add outflow from upstream to inflow
        sim_inflow .+= sn[i].outflow
    end

    # Run this node
    run_node!(
        sn[node_id], climate;
        inflow=sim_inflow, extraction=extraction, exchange=exchange
    )

    return nothing
end


"""
    run_node!(
        node::NetworkNode, climate::Climate;
        inflow=nothing, extraction=nothing, exchange=nothing
    )::Nothing

Run a specific node, and only that node, for all time steps.

# Arguments
- `node` : Any Streamfall NetworkNode
- `climate` : Climate data
- `inflow::Union{DataFrame, Vector, Number}` : Inflow to node
- `extraction::Union{DataFrame, Vector, Number}` : Extractions from this subcatchment
- `exchange::Union{DataFrame, Vector, Number}` : Groundwater flux
"""
function run_node!(
    node::NetworkNode, climate::Climate;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Nothing
    timesteps = sim_length(climate)
    prep_state!(node, timesteps)

    P_and_ET = climate_values(node, climate)

    node_name = node.name
    @inbounds for ts in 1:timesteps
        # Check if values are provided, otherwise default to 0.
        ext = timestep_value(ts, node_name, "extraction", extraction)
        flux = timestep_value(ts, node_name, "exchange", exchange)
        in_flow = timestep_value(ts, node_name, "inflow", inflow)

        run_timestep!(
            node, P_and_ET[ts, 1], P_and_ET[ts, 2], ts;
            inflow=in_flow, extraction=ext, exchange=flux
        )
    end

    return nothing
end

"""
    run_node!(
        node::NetworkNode, climate::Climate, ts::Int;
        inflow=nothing, extraction=nothing, exchange=nothing
    )

Run a specific node for a specified time step.

# Arguments
- `node` : A node in a network
- `climate` : Climate dataset
- `ts` : current time step
- `inflow` : Time series of inflows from upstream nodes.
- `extraction` : Time series of water orders (expects column of `_releases`)
- `exchange` : Time series of groundwater flux
"""
function run_node!(
    node::NetworkNode, climate::Climate, ts::Int;
    inflow=nothing, extraction=nothing, exchange=nothing
)
    node_name = node.name
    rain, et = climate_values(node, climate, ts)

    return run_timestep!(
        node, rain, et, ts;
        inflow=inflow, extraction=extraction, exchange=exchange
    )
end

include("Analysis/Analysis.jl")
include("plotting.jl")

# Nodes
export NetworkNode, GenericNode, GenericDirectNode
export IHACRES, IHACRESNode, IHACRESBilinearNode, ExpuhNode, DamNode, Climate
export create_node, GR4JNode, HyModNode, SimpleHyModNode, SYMHYDNode
export EnsembleNode, WeightedEnsembleNode

# Network
export find_inlets_and_outlets, inlets, outlets
export climate_values, node_names, get_node, get_node_id, get_prop, set_prop!
export param_info, update_params!, sim_length, reset!, prep_state!
export run_catchment!, run_basin!, run_node!, run_node_with_temp!
export run_timestep!
export calibrate!, calibrate_instances!, calibrate_weights!
export create_network, load_network, save_network

export best_candidate, best_fitness, best_params

# Data
export extract_flow, extract_climate, align_time_frame

# Plotting methods
export quickplot, plot_network, save_figure, temporal_cross_section

# Data interface (climate)
export timesteps

# Analysis
export Analysis

export apply_bias_correction

end  # end module
