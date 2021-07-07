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

# macro set_size(arr, len)


include("Network.jl")
include("Node.jl")
include("Climate.jl")
include("IHACRESNode.jl")
include("IHACRESExpuhNode.jl")
include("HyModNode.jl")
include("DamNode.jl")



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

    node = sn[node_id]
    if checkbounds(Bool, node.outflow, ts)
        if node.outflow[ts] != undef
            # already ran for this time step so no need to recurse further
            return node.outflow[ts], node.level[ts]
        end
    end
    
    sim_inflow = 0.0
    ins = inneighbors(sn.g, node_id)
    if !isempty(ins)
        for i in ins
            # Get inflow from previous node
            res = run_node!(sn, i, climate, ts;
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

    ts_inflow = timestep_value(ts, node.name, "inflow", inflow)
    ts_inflow += sim_inflow

    run_func! = get_prop(sn, node_id, :nfunc)

    # Run for a time step, dependent on previous state
    return run_func!(node, climate, ts; inflow=ts_inflow, extraction=extraction, exchange=exchange)
end


"""
    run_basin!(sn::StreamfallNetwork, climate::Climate; extraction=nothing, exchange=nothing)

Run scenario for an entire catchment/basin.
"""
function run_basin!(sn::StreamfallNetwork, climate::Climate;
                    inflow=nothing, extraction=nothing, exchange=nothing)
    _, outlets = find_inlets_and_outlets(sn)
    for outlet in outlets
        run_node!(sn, outlet, climate;
                  inflow=inflow, extraction=extraction, exchange=exchange)
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
function run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate; 
                   inflow=nothing, extraction=nothing, exchange=nothing)
    func = MetaGraphs.get_prop(sn.mg, node_id, :nfunc)
    func(sn, node_id, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    # timesteps = sim_length(climate)
    # for ts in 1:timesteps
    #     func(sn, node_id, climate, ts;
    #         inflow=inflow, extraction=extraction, exchange=exchange)
    # end
end


"""
    run_node!(node::NetworkNode, climate;
              inflow=nothing, extraction=nothing, exchange=nothing)

Run a specific node, and only that node, for all time steps.

# Arguments
- `node::NetworkNode` : Any Streamfall NetworkNode
- `climate::Climate` : Climate data
- `inflow::Union{DataFrame, Vector, Number}` : Inflow to node
- `extraction::Union{DataFrame, Vector, Number}` : Extractions from this subcatchment
- `exchange::Union{DataFrame, Vector, Number}` : Groundwater flux
"""
function run_node!(node::NetworkNode, climate; inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    # prep_state!(node, timesteps)

    node_name = node.name
    for ts in 1:timesteps
        rain, et = climate_values(node, climate, ts)
        ext = timestep_value(ts, node_name, "extraction", extraction)
        flux = timestep_value(ts, node_name, "exchange", exchange)
        in_flow = timestep_value(ts, node_name, "inflow", inflow)

        run_node!(node, rain, et, ts; inflow=in_flow, extraction=ext, exchange=flux)
    end

    return node.outflow, node.level
end


include("metrics.jl")
include("calibration.jl")


export @def
export NetworkNode, GenericNode, GenericDirectNode
export IHACRES, IHACRESNode, BilinearNode, ExpuhNode, DamNode, Climate
export HyModNode, SimpleHyModNode
export find_inlets_and_outlets, inlets, outlets, create_network, create_node
export climate_values, get_node, get_node_id, get_prop, set_prop!
export param_info, update_params!, sim_length, reset!
export run_catchment!, run_basin!, run_node!, run_node_with_temp!
export calibrate!, plot_network, save_figure

end  # end module
