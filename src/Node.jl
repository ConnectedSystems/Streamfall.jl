using DataFrames
using MetaGraphs


abstract type NetworkNode end

@def network_node begin
    node_id::String
    area::Float64  # area in km^2
    route::Bool
end


"""
Retrieve network node_id for a given gauge (by name).
"""
function get_node_id(mg::MetaGraph, gauge_id::String)::Int64
    v = collect(MetaGraphs.filter_vertices(mg, :name, gauge_id))
    @assert length(v) == 1 || "Found more than 1 node with gauge $(gauge_id)"
    return v[1]
end


"""
Retrieve node_id and node property for a specified gauge.

# Arguments
- sn::StreamfallNetwork
- gauge_id::String : name/id for node of interest

# Returns
Tuple, node position id and node object
"""
get_gauge(sn::StreamfallNetwork, gauge_id::String) = get_gauge(sn.mg, gauge_id)
function get_gauge(mg, gauge_id::String)::Tuple
    v_id = get_node_id(mg, gauge_id)
    return v_id, MetaGraphs.get_prop(mg, v_id, :node)
end


get_node(sn::StreamfallNetwork, v_id::Int64)::NetworkNode = MetaGraphs.get_prop(sn.mg, v_id, :node)
get_node_name(sn::StreamfallNetwork, v_id) = MetaGraphs.get_prop(sn.mg, v_id, :name)

get_node(mg, v_id) = MetaGraphs.get_prop(mg, v_id, :node)
get_node_name(mg, v_id) = MetaGraphs.get_prop(mg, v_id, :name)
