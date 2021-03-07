module Waterflow

using LightGraphs, MetaGraphs, Distributed, DataFrames
using Infiltrator


# Can't use string, DLL location has to be a const
# (which makes sense but still, many hours wasted!)
# https://github.com/JuliaLang/julia/issues/29602
const ihacres = "../../ihacres_nim/lib/ihacres.dll"

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


include("Node.jl")
include("StreamNode.jl")
include("DamNode.jl")
include("Climate.jl")

"""Run a model attached to a node.

water_order: Volume of water to be extracted (in ML/timestep)
exchange: Volume of flux (in ML/timestep), where negative values are losses to the groundwater system.
"""
function run_node!(mg::MetaGraph, g::AbstractGraph, node_id::Int, climate::Climate, timestep::Int; 
                   water_order=nothing, exchange=nothing)

    curr_node = get_prop(mg, node_id, :node)
    if checkbounds(Bool, curr_node.outflow, timestep)
        # already ran for this time step so no need to recurse further
        return curr_node.outflow[timestep], curr_node.level[timestep]
    end

    outflow = 0.0
    inflow = 0.0

    ins = inneighbors(g, node_id)
    if length(ins) == 0
        inflow = 0.0
        src_name = "Out-of-Catchment"
    else
        inflow = 0.0
        for i in ins
            src_name = get_prop(mg, i, :name)
            # Get inflow from previous node
            upstream_flow, upstream_level = run_node!(mg, g, i, climate, timestep)
            inflow += upstream_flow
        end
    end

    gauge_id = get_prop(mg, node_id, :name)
    rain, et = climate_values(curr_node, climate, timestep)
    if ismissing(rain) | ismissing(et)
        @warn "No climate data found for node $(node_id) ($(gauge_id)) using column markers $(climate.rainfall_id) and $(climate.et_id)"
        return 0.0
    end

    wo = 0.0
    if !isnothing(water_order)
        release_col = filter(x -> occursin(gauge_id, string(x))
                                  & occursin("releases", string(x)),
                                  names(water_order))
        
        wo = checkbounds(Bool, water_order.Date, timestep) ? water_order[timestep, release_col][1] : 0.0
    end

    ex = 0.0
    if !isnothing(exchange)
        exchange_col = filter(x -> occursin(gauge_id, x)
                                  & occursin("exchange", x),
                                  names(exchange))
        ex = checkbounds(Bool, exchange.Date, timestep) ? exchange[timestep, exchange_col][1] : 0.0
    end

    # Calculate outflow for this node
    func = get_prop(mg, node_id, :nfunc)
    if curr_node isa StreamNode
        outflow, level = func(curr_node, rain, et, inflow, 0.0)
    elseif curr_node isa DamNode
        outflow, level = func(curr_node, rain, et, inflow, wo, ex)
    else
        throw(ArgumentError("Unknown node type!"))
    end

    return outflow, level
end

include("Network.jl")


export @def
export ihacres, StreamNode, DamNode, Climate
export find_inlets_and_outlets, create_network, create_node
export update_params!, run_node!, reset!, sim_length
export climate_values

end  # end module
