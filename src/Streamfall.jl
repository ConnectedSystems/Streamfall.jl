module Streamfall

using Statistics
using LightGraphs, MetaGraphs, Distributed, DataFrames


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


"""@def macro

Inline code to avoid repetitious declarations.
"""
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


@def add_preprefix begin
    if !isnothing(id_prefix)
        prefix = id_prefix * prefix
    end
end


include("Network.jl")
include("Node.jl")
include("IHACRESNode.jl")
include("IHACRESExpuhNode.jl")
include("HyModNode.jl")
include("DamNode.jl")
include("Climate.jl")


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
- `extraction::DataFrame` : Volume of water to be extracted (in ML/timestep)
- `exchange::DataFrame` : Volume of flux (in ML/timestep), where negative values are losses to the groundwater system
"""
function run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate, timestep::Int;
                   extraction::Union{DataFrame, Nothing}=nothing,
                   exchange::Union{DataFrame, Nothing}=nothing)
    mg, g = sn.mg, sn.g

    curr_node = MetaGraphs.get_prop(mg, node_id, :node)
    if checkbounds(Bool, curr_node.outflow, timestep)
        if curr_node.outflow[timestep] != undef
            # already ran for this time step so no need to recurse further
            return curr_node.outflow[timestep], curr_node.level[timestep]
        end
    end

    outflow = 0.0
    inflow = 0.0

    ins = inneighbors(g, node_id)
    if isempty(ins)
        inflow = 0.0
    else
        inflow = 0.0
        for i in ins
            # Get inflow from previous node
            upstream_flow, _ = run_node!(sn, i, climate, timestep; extraction=extraction, exchange=exchange)
            inflow += upstream_flow
        end
    end

    gauge_id = curr_node.node_id
    rain, et = climate_values(curr_node, climate, timestep)
    wo = timestep_value(timestep, gauge_id, "releases", extraction)
    ex = timestep_value(timestep, gauge_id, "exchange", exchange)

    # Calculate outflow for this node
    func = MetaGraphs.get_prop(mg, node_id, :nfunc)
    outflow, level = func(curr_node, rain, et, inflow, wo, ex)

    return outflow, level
end


"""
    run_basin!(sn::StreamfallNetwork, climate::Climate; extraction=nothing, exchange=nothing)

Run scenario for an entire catchment/basin.
"""
function run_basin!(sn::StreamfallNetwork, climate::Climate;
                    extraction=nothing, exchange=nothing)
    _, outlets = find_inlets_and_outlets(sn)
    for outlet in outlets
        run_node!(sn, outlet, climate; extraction=extraction, exchange=exchange)
    end
end

run_catchment! = run_basin!


"""
    run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate;
              extraction=nothing, exchange=nothing)::Nothing

Run model for all time steps, recursing upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id::Int` : node to run in the network
- `climate::Climate` : Climate object holding rainfall and evaporation data (or temperature)
- `extraction::DataFrame` : water orders for each time step (defaults to nothing)
- `exchange::DataFrame` : exchange with groundwater system at each time step (defaults to nothing)
"""
function run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate; extraction=nothing, exchange=nothing)::Nothing
    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(sn, node_id, climate, ts; extraction=extraction, exchange=exchange)
    end
end


"""
    run_node!(node::NetworkNode, climate;
              inflow=nothing, extraction=nothing, exchange=nothing)

Run a specific node, and only that node, for all time steps.

# Arguments
- `node::NetworkNode` :
- `climate::Climate` :
- `inflow::Vector` : Time series of inflows.
- `extraction::Vector` : Time series of water extractions from this subcatchment
- `exchange::Vector` : Time series of groundwater flux
"""
function run_node!(node::NetworkNode, climate; inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    gauge_id = node.node_id

    for ts in (1:timesteps)
        run_node!(node, climate, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow, node.level
end


"""
    run_node!(node::NetworkNode, climate::Climate, ts::Int; 
              inflow=nothing, extraction=nothing, exchange=nothing)

Run a specific node for a specified time step.

# Arguments
- `node::NetworkNode` :
- `climate::Climate` :
- `inflow::DataFrame` : Time series of inflows from any upstream node.
- `extraction::DataFrame` : Time series of water extractions from this subcatchment
- `exchange::DataFrame` : Time series of groundwater flux
"""
function run_node!(node::NetworkNode, climate::Climate, ts::Int; inflow=nothing, extraction=nothing, exchange=nothing)

    gauge_id = node.node_id
    if checkbounds(Bool, node.outflow, ts)
        if curr_node.outflow[ts] != undef
            # already ran for this time step so no need to recurse further
            return curr_node.outflow[ts], curr_node.level[ts]
        end
    end

    rain, et = climate_values(node, climate, ts)
    wo = timestep_value(ts, gauge_id, "extraction", extraction)
    ex = timestep_value(ts, gauge_id, "exchange", exchange)
    in_flow = timestep_value(ts, gauge_id, "inflow", inflow)

    run_node!(node, rain, et, in_flow, wo, ex)

    return node.outflow[ts], node.level[ts]
end



include("metrics.jl")


export @def
export IHACRES, IHACRESNode, BilinearNode, ExpuhNode, DamNode, Climate
export HyModNode, SimpleHyModNode
export find_inlets_and_outlets, inlets, outlets, create_network, create_node
export param_info, update_params!, run_node!, reset!, sim_length, run_catchment!
export run_node_with_temp!
export climate_values, get_node, get_gauge, get_node_id, get_prop
export set_prop!

end  # end module
