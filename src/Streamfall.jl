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
include("DamNode.jl")
include("Climate.jl")


function timestep_value(timestep::Int, gauge_id::String, col_partial::String, dataset=nothing)::Float64
    amount::Float64 = 0.0
    if !isnothing(dataset)
        target_col = filter(x -> occursin(gauge_id, x)
                                & occursin(col_partial, x),
                                names(dataset))
        if !isempty(target_col)
            amount = checkbounds(Bool, dataset.Date, timestep) ? dataset[timestep, target_col][1] : 0.0
        end
    end

    return amount
end


"""Find common time frame between time series."""
function find_common_timeframe(timeseries::T...) where {T<:DataFrame}
    for t in timeseries
        @assert "Date" in names(t)
    end

    min_date = maximum([x.Date[1] for x in timeseries])
    max_date = minimum([x.Date[end] for x in timeseries])
    return (min_date, max_date)
end


"""Run a model attached to a node for a given time step.
Recurses upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id::Int`
- `climate::Climate`
- `timestep::Int`
- `water_order::Vector` : Volume of water to be extracted (in ML/timestep)
- `exchange::Vector` : Volume of flux (in ML/timestep), where negative values are losses to the groundwater system
"""
run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate, timestep::Int; water_order=nothing, exchange=nothing) =
    run_node!(sn.mg, sn.g, node_id, climate, timestep; water_order=water_order, exchange=exchange)
function run_node!(mg::MetaGraph, g::AbstractGraph, node_id::Int, climate::Climate, timestep::Int; 
                   water_order=nothing, exchange=nothing)

    curr_node = MetaGraphs.get_prop(mg, node_id, :node)
    if checkbounds(Bool, curr_node.outflow, timestep)
        # already ran for this time step so no need to recurse further
        return curr_node.outflow[timestep], curr_node.level[timestep]
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
            upstream_flow, upstream_level = run_node!(mg, g, i, climate, timestep; water_order=water_order, exchange=exchange)
            inflow += upstream_flow
        end
    end

    gauge_id = curr_node.node_id
    rain, et = climate_values(curr_node, climate, timestep)
    if ismissing(rain) | ismissing(et)
        @warn "No climate data found for node $(node_id) ($(gauge_id)) using column markers $(climate.rainfall_id) and $(climate.et_id)"
        return 0.0
    end

    wo = timestep_value(timestep, gauge_id, "releases", water_order)
    ex = timestep_value(timestep, gauge_id, "exchange", exchange)

    # Calculate outflow for this node
    func = MetaGraphs.get_prop(mg, node_id, :nfunc)
    if curr_node isa IHACRESNode
        outflow, level = func(curr_node, rain, et, inflow, wo, ex)
    elseif curr_node isa ExpuhNode
        outflow, level = func(curr_node, rain, et, inflow, wo, ex)
    elseif curr_node isa DamNode
        outflow, level = func(curr_node, rain, et, inflow, wo, ex)
    else
        throw(ArgumentError("Unknown node type!"))
    end

    return outflow, level
end


"""Run scenario for an entire catchment."""
function run_catchment!(sn::StreamfallNetwork, climate; water_order=nothing, exchange=nothing) 
    run_catchment!(sn.mg, sn.g, climate; water_order=water_order, exchange=exchange)
end
function run_catchment!(mg, g, climate; water_order=nothing, exchange=nothing)
    inlets, outlets = find_inlets_and_outlets(g)
    for outlet in outlets
        run_node!(mg, g, outlet, climate; water_order=water_order, exchange=exchange)
    end
end


"""Run model for all time steps, recursing upstream as needed.

# Arguments
- `sn::StreamfallNetwork`
- `node_id::Int` : node to run in the network
- `climate::Climate` : Climate object holding rainfall and evaporation data (or temperature)
- `water_order::Vector` : water orders for each time step (defaults to nothing)
- `exchange::Vector` : exchange with groundwater system at each time step (defaults to nothing)
"""
function run_node!(sn::StreamfallNetwork, node_id, climate; water_order=nothing, exchange=nothing)::Nothing
    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(sn.mg, sn.g, node_id, climate, ts; water_order=water_order, exchange=exchange)
    end
end


function run_node!(mg, g, target_node, climate; water_order=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(mg, g, target_node, climate, ts; water_order=water_order, exchange=exchange)
    end
end


"""
Run a specific node, and only that node, for all time steps.

# Arguments
- `node::NetworkNode` :
- `climate::Climate` :
- `inflow::Vector` : Time series of inflows from any upstream node.
- `water_order::Vector` : Time series of water extractions from this subcatchment
- `exchange::Vector` : Time series of groundwater flux
"""
function run_node!(node::NetworkNode, climate; inflow=nothing, water_order=nothing, exchange=nothing)
    timesteps = sim_length(climate)

    for ts in (1:timesteps)
        if checkbounds(Bool, curr_node.outflow, ts)
            # already ran for this time step so no need to recurse further
            return curr_node.outflow[ts], curr_node.level[ts]
        end

        gauge_id = node.node_id
        rain, et = climate_values(node, climate, ts)
        if ismissing(rain) | ismissing(et)
            @warn "No climate data found for node $(node_id) ($(gauge_id)) using column markers $(climate.rainfall_id) and $(climate.et_id)"
            return 0.0
        end

        wo = 0.0
        ex = 0.0
        if !isnothing(water_order)
            wo = water_order[ts]
        end

        if !isnothing(exchange)
            ex = exchange[ts]
        end

        if !isnothing(inflow)
            if inflow isa Vector
                in_flow = inflow[ts]
            elseif inflow isa Number
                in_flow = inflow
            end
        end

        run_node!(node, rain, et, in_flow, wo, ex)
    end

    return node.outflow, node.level
end


include("metrics.jl")


export @def
export IHACRES, IHACRESNode, BilinearNode, ExpuhNode, DamNode, Climate
export find_inlets_and_outlets, inlets, outlets, create_network, create_node
export param_info, update_params!, run_node!, reset!, sim_length, run_catchment!
export run_node_with_temp!
export climate_values, get_node, get_gauge, get_node_id, get_prop
export set_prop!

end  # end module
