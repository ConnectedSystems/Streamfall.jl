using DataFrames
using PrettyTables
import MetaGraphs


abstract type NetworkNode end

@def network_node begin
    name::String
    area::Float64
end


"""
    param_info(node::NetworkNode)

Generic parameter information extractor.

Extracts parameter names, values, and bounds
"""
function param_info(node::NetworkNode; kwargs...)::Tuple
    tmp = Model(node)
    values = collect(tmp.val)
    bounds = collect(tmp.bounds)
    param_names = collect(tmp.fieldname)
    
    return param_names, values, bounds
end


"""Retrieve network node_id for a given gauge (by name)."""
function get_node_id(mg::MetaGraph, node_name::String)::Int64
    v = collect(MetaGraphs.filter_vertices(mg, :name, node_name))
    @assert length(v) > 0 || error("Node with name '$(node_name)' not found")
    @assert length(v) == 1 || error("Found more than 1 node with name '$(node_name)'")
    return v[1]
end


"""
    get_node(sn::StreamfallNetwork, node_name::String)

Retrieve node_id and node property for a specified gauge.

# Arguments
- sn::StreamfallNetwork
- node_name::String : name of node of interest

# Returns
Tuple, node position id and node object
"""
get_node(sn::StreamfallNetwork, node_name::String) = get_node(sn.mg, node_name)
function get_node(mg, node_name::String)::Tuple
    v_id = get_node_id(mg, node_name)
    return v_id, MetaGraphs.get_prop(mg, v_id, :node)
end


get_node(sn::StreamfallNetwork, v_id::Int)::NetworkNode = MetaGraphs.get_prop(sn.mg, v_id, :node)
get_node_name(sn::StreamfallNetwork, v_id::Int) = MetaGraphs.get_prop(sn.mg, v_id, :name)

get_node(mg, v_id::Int) = MetaGraphs.get_prop(mg, v_id, :node)
get_node_name(mg, v_id::Int) = MetaGraphs.get_prop(mg, v_id, :name)


"""
    subcatchment_data(node::NetworkNode, data::DataFrame)

Extract data for a given node from a dataframe.
"""
function subcatchment_data(node::NetworkNode, data::DataFrame, partial::String=nothing)::DataFrame
    name = node.name
    cols = filter(x -> occursin(name, string(x)), names(data))

    if !isnothing(partial)
        cols = filter(x -> occursin(name, x), cols)
    end

    return data[:, cols]
end


function extract_node_spec(node::NetworkNode)
    spec = Dict()

    area = node.area
    param_names, x0, _ = param_info(node)
    params = Dict(zip(param_names, x0))

    node_type = typeof(node).name
    spec = Dict(
        "node_type" => node_type,
        "area" => area,
        "parameters" => params
    )

    return spec
end


Base.show(io::IO, ::MIME"text/plain", n::NetworkNode) = show(io, n)
function Base.show(io::IO, n::NetworkNode)
    println(io, "Name: $(n.name) [$(typeof(n).name)]")
    println(io, "Area: $(n.area)")

    param_names, x0, bounds = param_info(n; with_level=false)
    lb, ub = zip(bounds...)
    details = hcat(param_names, x0, [lb...], [ub...])

    pretty_table(io, details, ["Parameter", "Value", "Lower Bound", "Upper Bound"])
end
