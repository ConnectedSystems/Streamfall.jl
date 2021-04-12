module Streamfall

using LightGraphs, MetaGraphs, Distributed, DataFrames
using Statistics


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


include("Node.jl")
include("IHACRESNode.jl")
include("DamNode.jl")
include("Climate.jl")


function timestep_value(timestep, gauge_id::String, col_partial::String, dataset=nothing)
    amount = 0.0
    if !isnothing(dataset)
        target_col = filter(x -> occursin(gauge_id, string(x))
                                & occursin(col_partial, string(x)),
                                names(dataset))
        if !isempty(target_col)
            amount = checkbounds(Bool, dataset.Date, timestep) ? dataset[timestep, target_col][1] : 0.0
        end
    end

    return amount
end

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
    if isempty(ins)
        inflow = 0.0
        src_name = "Out-of-Catchment"
    else
        inflow = 0.0
        for i in ins
            src_name = get_prop(mg, i, :name)
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
    func = get_prop(mg, node_id, :nfunc)
    if curr_node isa IHACRESNode
        outflow, level = func(curr_node, rain, et, inflow, 0.0)
    elseif curr_node isa DamNode
        outflow, level = func(curr_node, rain, et, inflow, wo, ex)
    else
        throw(ArgumentError("Unknown node type!"))
    end

    return outflow, level
end


function run_catchment!(mg, g, climate; water_order=nothing, exchange=nothing)
    inlets, outlets = find_inlets_and_outlets(g)
    for outlet in outlets
        run_node!(mg, g, outlet, climate; water_order=water_order, exchange=exchange)
    end
end


function run_node!(mg, g, target_node, climate; water_order=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(mg, g, target_node, climate, ts; water_order=water_order, exchange=exchange)
    end
end


include("Network.jl")
include("metrics.jl")


export @def
export IHACRES, IHACRESNode, DamNode, Climate
export find_inlets_and_outlets, create_network, create_node
export param_info, update_params!, run_node!, reset!, sim_length, run_catchment!
export climate_values, get_node, get_gauge, get_node_id

end  # end module
