using LightGraphs, MetaGraphs
using ModelParameters
# using Streamfall


struct StreamfallNetwork
    mg::MetaGraph
    g::SimpleDiGraph
end


function set_prop!(sn::StreamfallNetwork, nid::Int64, prop::Symbol, var::Any)::Nothing
    MetaGraphs.set_prop!(sn.mg, nid, prop, var)

    return nothing
end


function get_prop(sn::StreamfallNetwork, nid::Int64, prop::Symbol)::Any
    return MetaGraphs.get_prop(sn.mg, nid, prop)
end


"""Determine a node's connection"""
function in_or_out(G, v)
    ins = length(inneighbors(G, v))
    outs = length(outneighbors(G, v))

    inlet = false
    outlet = false
    if ins == 0
        inlet = true
    elseif outs == 0
        outlet = true
    end

    return v, inlet, outlet
end


"""Find all inlets and outlets in a network."""
find_inlets_and_outlets(sn::StreamfallNetwork) = find_inlets_and_outlets(sn.g)
function find_inlets_and_outlets(G)
    vs = vertices(G)
    num_vs::Int64 = length(vs)

    ins_outs = pmap(in_or_out, repeat([G], num_vs), vs)
    inlets = Int64[]
    outlets = Int64[]
    for row in ins_outs
        if (row[2] == false) & (row[3] == false)
            continue
        end

        if row[2] == true
            push!(inlets, row[1])
        elseif row[3] == true
            push!(outlets, row[1])
        end
    end

    return inlets, outlets
end


"""Create a node if needed"""
create_node(sn::StreamfallNetwork, node_id, details, nid) = create_node(sn.mg, node_id, details, nid)
function create_node(mg, node_id, details, nid)
    details = copy(details)
    match = collect(MetaGraphs.filter_vertices(mg, :name, node_id))
    if isempty(match)
        node_type = details["node_type"]

        dtype = eval(Symbol(node_type))
        n = nothing
        try
            n = dtype(node_id, details)
        catch
            throw(ArgumentError("Unsupported node type: $(node_type)"))
        end

        set_props!(mg, nid, Dict(:name=>node_id,
                                 :node=>n,
                                 :nfunc=>run_node!))
        
        this_id = nid
        nid = nid + 1
    else
        this_id = match[1]
    end

    return this_id, nid
end


"""Create a network from a YAML-derived spec

    network = YAML.load_file("example_network.yml")
    sn = create_network("Example Network", network)

"""
function create_network(name::String, network::Dict)
    num_nodes = length(network)
    g = SimpleDiGraph(num_nodes)
    mg = MetaGraph(g)
    MetaGraphs.set_prop!(mg, :description, name)
    
    nid = 1
    for (node, details) in network
        node_id = string(node)

        this_id, nid = create_node(mg, node_id, details, nid)

        inlets = details["inlets"]
        in_id = nid
        if !isnothing(inlets)
            for inlet in inlets
                in_id, nid = create_node(mg, string(inlet), network[inlet], nid)
                add_edge!(g, in_id, this_id)
            end
        end

        outlets = details["outlets"]
        out_id = in_id
        if !isnothing(outlets)
            for outlet in outlets
                out_id, nid = create_node(mg, string(outlet), network[outlet], nid)
                add_edge!(g, this_id, out_id)
            end
        end
    end

    sn = StreamfallNetwork(mg, g)

    return sn
end


"""
Reset a network.
"""
reset!(sn::StreamfallNetwork) = reset!(sn.mg, sn.g)
function reset!(mg, g)::Nothing
    v_ids = vertices(g)
    for i in v_ids
        curr_node = MetaGraphs.get_prop(mg, i, :node)
        reset!(curr_node)
    end
end