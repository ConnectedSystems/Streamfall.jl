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
include("Nodes/IHACRES/IHACRESNode.jl")
include("Nodes/IHACRES/IHACRESExpuhNode.jl")
include("Nodes/GR4J/GR4JNode.jl")
include("Nodes/HyMod/HyModNode.jl")
include("Nodes/SYMHYD/SYMHYDNode.jl")
include("Nodes/DamNode.jl")
include("Nodes/EnsembleNode.jl")


function timestep_value(timestep::Int, gauge_id::String, col_partial::String,
                        dataset::Union{DataFrame, Vector, Number, Nothing}=nothing)::Float64
    amount::Float64 = 0.0
    if !isnothing(dataset)
        if dataset isa Vector
            amount = dataset[ts]
        elseif dataset isa DataFrame
            # Extract data for time step based on partial match
            target_col = filter(x -> occursin(gauge_id, x)
                                & occursin(col_partial, x),
                                names(dataset))
            if !isempty(target_col)
                amount = checkbounds(Bool, dataset.Date, timestep) ? dataset[timestep, target_col][1] : 0.0
            end
        elseif dataset isa Number
            amount = dataset
        end
    end

    return amount
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

    return modded
end


"""
    run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate, timestep::Int;
              extraction::Union{DataFrame, Nothing}=nothing,
              exchange::Union{DataFrame, Nothing}=nothing)

Run a model attached to a node for a given time step.
Recurses upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id::Int`
- `climate::Climate`
- `timestep::Int`
- `inflow::DataFrame` : Additional inflow to consider (in ML/timestep)
- `extraction::DataFrame` : Volume of water to be extracted (in ML/timestep)
- `exchange::DataFrame` : Volume of flux (in ML/timestep), where negative values are losses to the groundwater system
"""
function run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate, timestep::Int;
                   inflow::Union{DataFrame, Nothing}=nothing,
                   extraction::Union{DataFrame, Nothing}=nothing,
                   exchange::Union{DataFrame, Nothing}=nothing)
    ts = timestep

    sim_inflow = 0.0
    ins = inneighbors(sn.mg, node_id)

    if !isempty(ins)
        for i in ins
            # Get inflow from previous node
            res = run_node!(sn, i, climate, timestep;
                            inflow=inflow,
                            extraction=extraction,
                            exchange=exchange)
            if res isa Number
                sim_inflow += res
            elseif length(res) > 1
                # get outflow from (outflow, level)
                sim_inflow += res[1]
            end
        end
    end

    node = sn[node_id]
    ts_inflow = timestep_value(ts, node.name, "inflow", inflow)
    ts_inflow += sim_inflow

    run_func! = get_prop(sn, node_id, :nfunc)

    # Run for a time step, dependent on previous state
    run_func!(node, climate, ts; inflow=ts_inflow, extraction=extraction, exchange=exchange)
end


"""
    run_basin!(sn::StreamfallNetwork, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)

Run scenario for an entire catchment/basin.
"""
function run_basin!(sn::StreamfallNetwork, climate::Climate;
                    inflow=nothing, extraction=nothing, exchange=nothing)
    _, outlets = find_inlets_and_outlets(sn)
    @inbounds for outlet in outlets
        run_node!(sn, outlet, climate;
                  inflow=inflow, extraction=extraction, exchange=exchange)
    end
end

run_catchment! = run_basin!


"""
    run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate;
              extraction=nothing, exchange=nothing)::Nothing

Generic run method that runs a model attached to a given node for all time steps.
Recurses upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id::Int` : node to run in the network
- `climate::Climate` : Climate object holding rainfall and evaporation data (or temperature)
- `extraction::DataFrame` : water orders for each time step (defaults to nothing)
- `exchange::DataFrame` : exchange with groundwater system at each time step (defaults to nothing)
"""
function run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate;
                   inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    @inbounds for ts in 1:timesteps
        run_node!(sn, node_id, climate, ts;
                  inflow=inflow, extraction=extraction, exchange=exchange)
    end
end


"""
    run_node!(node::NetworkNode, climate::Climate;
              inflow=nothing, extraction=nothing, exchange=nothing)

Run a specific node, and only that node, for all time steps.

# Arguments
- `node::NetworkNode` : Any Streamfall NetworkNode
- `climate` : Climate data
- `inflow::Union{DataFrame, Vector, Number}` : Inflow to node
- `extraction::Union{DataFrame, Vector, Number}` : Extractions from this subcatchment
- `exchange::Union{DataFrame, Vector, Number}` : Groundwater flux
"""
function run_node!(node::NetworkNode, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)::Nothing
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


include("metrics.jl")
include("calibration.jl")

include("Analysis/Analysis.jl")
include("plotting.jl")

# Nodes
export NetworkNode, GenericNode, GenericDirectNode
export IHACRES, IHACRESNode, BilinearNode, ExpuhNode, DamNode, Climate
export create_node, GR4JNode, HyModNode, SimpleHyModNode, SYMHYDNode
export EnsembleNode, BaseEnsemble
export run_step!, run_timestep!

# Network
export find_inlets_and_outlets, inlets, outlets, create_network, create_node
export climate_values, get_node, get_node_id, get_prop, set_prop!
export param_info, update_params!, sim_length, reset!
export run_catchment!, run_basin!, run_node!, run_node_with_temp!
export calibrate!

# Data
export extract_flow, extract_climate

# plotting methods
export quickplot, plot_network, save_figure

# data interface (climate)
export timesteps

# Analysis
export Analysis

end  # end module
