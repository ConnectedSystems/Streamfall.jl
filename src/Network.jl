using LightGraphs, MetaGraphs
using ModelParameters
using Streamfall


struct StreamfallNetwork
    mg::MetaGraph
    g::SimpleDiGraph
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
function create_node(mg, node_id, details, nid)
    details = copy(details)
    match = collect(filter_vertices(mg, :name, node_id))
    if isempty(match)
        node_type = details["node_type"]

        # TODO: Clean this up...
        # Needs to just dispatch on node_type and specification
        if node_type == "IHACRESNode"
            n = IHACRESNode(node_id, details)
        elseif node_type == "DamNode"
            n = DamNode(node_id, details)
        elseif node_type == "ExpuhNode"
            n = ExpuhNode(node_id, details)
        else
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
    mg, g = create_network("Example Network", network)

"""
function create_network(name::String, network::Dict)
    num_nodes = length(network)
    g = SimpleDiGraph(num_nodes)
    mg = MetaGraph(g)
    set_prop!(mg, :description, name)
    
    # sn = StreamfallNetwork(mg, g)
    
    nid = 1
    for (node, details) in network
        node_id = string(node)

        this_id, nid = create_node(mg, node_id, details, nid)
        # @info "Creating $(node_id) as $(this_id)"

        inlets = details["inlets"]
        in_id = nid
        if !isnothing(inlets)
            for inlet in inlets
                in_id, nid = create_node(mg, string(inlet), network[inlet], nid)
                # @info "Inlet: Creating link from $(inlet) to $(node_id) | $(in_id) -> $(this_id)"
                add_edge!(g, in_id, this_id)
            end
        end

        outlets = details["outlets"]
        out_id = in_id
        if !isnothing(outlets)
            for outlet in outlets
                out_id, nid = create_node(mg, string(outlet), network[outlet], nid)
                # @info "Outlet: Creating link from $(node_id) to $(outlet) | $(this_id) -> $(out_id)"
                add_edge!(g, this_id, out_id)
            end
        end
    end

    return mg, g
end


"""
Reset a network.
"""
function reset!(mg, g)
    v_ids = vertices(g)
    for i in v_ids
        curr_node = get_prop(mg, i, :node)
        reset!(curr_node)
    end
end