using Cairo, Compose
using LightGraphs, MetaGraphs, GraphPlot
using ModelParameters

import YAML: write_file


struct StreamfallNetwork
    mg::MetaGraph
    g::SimpleDiGraph
end


Base.getindex(sn::StreamfallNetwork, n::String) = get_node(sn, n)
Base.getindex(sn::StreamfallNetwork, nid::Int) = get_node(sn, nid)


function set_prop!(sn::StreamfallNetwork, nid::Int64, prop::Symbol, var::Any)::Nothing
    MetaGraphs.set_prop!(sn.mg, nid, prop, var)

    return nothing
end


function get_prop(sn::StreamfallNetwork, nid::Int64, prop::Symbol)::Any
    return MetaGraphs.get_prop(sn.mg, nid, prop)
end


"""Determine a node's connection"""
function in_or_out(g, v)
    ins = length(inneighbors(g, v))
    outs = length(outneighbors(g, v))

    inlet = false
    outlet = false
    if outs == 0
        outlet = true
    elseif ins == 0
        inlet = true
    end

    return v, inlet, outlet
end


"""Find all inlets and outlets in a network."""
function find_inlets_and_outlets(sn::StreamfallNetwork)
    g = sn.g
    vs = vertices(g)
    num_vs::Int64 = length(vs)

    ins_outs = pmap(in_or_out, repeat([g], num_vs), vs)
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


inlets(sn::StreamfallNetwork, nid::Number) = inneighbors(sn.g, nid)
outlets(sn::StreamfallNetwork, nid::Number) = outneighbors(sn.g, nid)


"""
    area(sn::StreamfallNetwork)::Float64

Total area represented by a network.
"""
function area(sn::StreamfallNetwork)::Float64
    num_nodes = nv(sn.g)
    area = 0.0
    for nid in 1:num_nodes
        area += sn[nid].area
    end

    return area
end


"""
    inlets(sn::StreamfallNetwork, node_name::String)

Find nodes which provides inflows for given node.
"""
function inlets(sn::StreamfallNetwork, node_name::String)::Array{Int}
    nid, _ = sn[node_name]
    return inneighbors(sn.g, nid)
end


"""
    outlets(sn::StreamfallNetwork, node_name::String)

Find node immediately downstream from given node.
"""
function outlets(sn::StreamfallNetwork, node_name::String)::Array{Int}
    nid, _ = sn[node_name]
    return outneighbors(sn.g, nid)
end


"""
    create_node(mg::MetaGraph, node_name::String, details::Dict, nid::Int)

Create a node specified with given name (if it does not exist).
"""
function create_node(mg::MetaGraph, node_name::String, details::Dict, nid::Int)
    details = copy(details)

    match = collect(MetaGraphs.filter_vertices(mg, :name, node_name))
    if isempty(match)
        node_type = details["node_type"]

        dtype = eval(Symbol(node_type))
        n = nothing
        try
            n = dtype(node_name, details)
        catch
            throw(ArgumentError("Unsupported node type: $(node_type)"))
        end

        # Set function for node if specified
        if haskey(details, "func")
            func_spec = details["func"]
            if func_spec isa Function
                func = func_spec
            elseif func_spec isa String
                # func = eval(Symbol(func_spec))
                func = eval(Meta.parse(func_spec))
            end
        else
            func = run_node!
        end

        set_props!(mg, nid, Dict(:name=>node_name,
                                :node=>n,
                                :nfunc=>func))
        
        this_id = nid
        nid = nid + 1
    else
        this_id = match[1]
    end

    return this_id, nid
end


"""
    create_network(name::String, network::Dict)::StreamfallNetwork

Create a StreamNetwork from a YAML-derived specification.

# Example
```julia-repl
julia> network_spec = YAML.load_file("example_network.yml")
julia> sn = create_network("Example Network", network_spec)
```
"""
function create_network(name::String, network::Dict)::StreamfallNetwork
    num_nodes = length(network)
    g = SimpleDiGraph(num_nodes)
    mg = MetaGraph(g)
    MetaGraphs.set_prop!(mg, :description, name)
    
    nid = 1
    for (node, details) in network
        n_name = string(node)

        this_id, nid = create_node(mg, n_name, details, nid)

        if haskey(details, "inlets")
            inlets = details["inlets"]
            in_id = nid
            if !isnothing(inlets)
                for inlet in inlets
                    in_id, nid = create_node(mg, string(inlet), network[inlet], nid)
                    add_edge!(g, in_id, this_id)
                end
            end
        end

        if haskey(details, "outlets")
            outlets = details["outlets"]
            out_id = in_id
            if !isnothing(outlets)
                msg = "Streamfall currently only supports a single outlet. ($(length(outlets)))"
                @assert length(outlets) <= 1 || throw(ArgumentError(msg))

                for outlet in outlets
                    out_id, nid = create_node(mg, string(outlet), network[outlet], nid)
                    add_edge!(g, this_id, out_id)
                end
            end
        end
    end

    sn = StreamfallNetwork(mg, g)

    return sn
end


"""
    reset!(sn::StreamfallNetwork)::Nothing

Reset a network.
"""
function reset!(sn::StreamfallNetwork)::Nothing
    mg, g = sn.mg, sn.g
    v_ids = vertices(g)
    for i in v_ids
        curr_node = MetaGraphs.get_prop(mg, i, :node)
        reset!(curr_node)
    end
end


function extract_node_spec!(sn::StreamfallNetwork, nid::Int, spec::Dict)::Nothing
    node = sn[nid]

    node_name = string(node.name)
    if haskey(spec, node_name)
        # This node already extracted.
        return
    end

    ins = inlets(sn, nid)
    outs = outlets(sn, nid)
    in_ids::Union{Array{String}, Nothing} = [sn[i].name for i in ins]
    out_ids::Union{Array{String}, Nothing} = [sn[i].name for i in outs]

    if length(out_ids) == 0
        out_ids = nothing
    end

    if length(in_ids) == 0
        in_ids = nothing
    else
        # Recurse upstream
        for n in ins
            extract_node_spec!(sn, n, spec)
        end
    end

    node_spec = extract_node_spec(node)
    network_spec = Dict(
        "inlets" => in_ids,
        "outlets" => out_ids
    )

    spec[node_name] = merge(node_spec, network_spec)

    return nothing
end


"""
    extract_network_spec(sn::StreamfallNetwork)

Extract network details
"""
function extract_network_spec(sn::StreamfallNetwork)::Dict
    _, outlets = find_inlets_and_outlets(sn)
    spec = Dict()
    for nid in outlets
        extract_node_spec!(sn, nid, spec)
    end

    # TODO: Make spec nicely ordered

    return spec
end


function save_network_spec(sn::StreamfallNetwork, fn::String)
    spec = extract_network_spec(sn)
    write_file(fn, spec)
end


Base.show(io::IO, ::MIME"text/plain", sn::StreamfallNetwork) = show(io, sn)
function Base.show(io::IO, sn::StreamfallNetwork)

    name = MetaGraphs.get_prop(sn.mg, :description)

    println(io, "Network Name: $(name)")
    println(io, "Represented Area: $(area(sn))")
    println("-"^17, "\n")

    vs = vertices(sn.g)

    for nid in vs
        println(io, "Node $(nid) : \n")
        show(io, sn[nid])
        print("\n")
    end
end


"""
    plot_network(sn::StreamfallNetwork)

Simple plot of stream network.
"""
function plot_network(sn::StreamfallNetwork; as_html=false)
    node_names = ["$(n.name)" for n in sn]

    if as_html
        plot_func = gplothtml
    else
        plot_func = gplot
    end

    plot_func(sn.g, nodelabel=node_names)
end


"""
    save_figure(sn::StreamfallNetwork, fn::String)

Save a figure of the network in SVG format.
"""
function save_figure(sn::StreamfallNetwork, fn::String)
    draw(SVG(fn, 16cm, 16cm), plot_network(sn))
end


function Base.iterate(sn::StreamfallNetwork, state=1)
    if state >= length(sn)
        return nothing
    end

    return (sn[state], state+1)
end


Base.length(sn::StreamfallNetwork) = nv(sn.g)



