using OrderedCollections

using Cairo, Compose
using Graphs, MetaGraphs, GraphPlot
using ModelParameters

import YAML: load_file, write_file


struct StreamfallNetwork
    mg::MetaDiGraph
end

function load_network(name::String, fn::String)
    network = load_file(
        fn;
        dicttype=OrderedDict
    )

    return create_network(name, network)
end


Base.getindex(sn::StreamfallNetwork, n::String) = get_node(sn, n)
Base.getindex(sn::StreamfallNetwork, nid::Int) = get_node(sn, nid)

function node_names(sn::StreamfallNetwork)
    verts = vertices(sn.mg)
    return [sn[v].name for v in verts]
end

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
function find_inlets_and_outlets(sn::StreamfallNetwork)::Tuple
    g = sn.mg
    vs = vertices(g)
    num_vs::Int64 = length(vs)

    ins_outs = map(in_or_out, repeat([g], num_vs), vs)
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


inlets(sn::StreamfallNetwork, nid::Number) = inneighbors(sn.mg, nid)
outlets(sn::StreamfallNetwork, nid::Number) = outneighbors(sn.mg, nid)


"""
    area(sn::StreamfallNetwork)::Float64

Total area represented by a network.
"""
function area(sn::StreamfallNetwork)::Float64
    num_nodes = nv(sn.mg)
    area = 0.0
    for nid in 1:num_nodes
        area += sn[nid].area
    end

    return area
end


"""
    inlets(sn::StreamfallNetwork, node_name::String)

Find ID(s) of nodes which provides inflows for given node.
"""
function inlets(sn::StreamfallNetwork, node_name::String)::Array{Int}
    nid, _ = sn[node_name]
    return inneighbors(sn.mg, nid)
end


"""
    outlets(sn::StreamfallNetwork, node_name::String)

Find ID(s) of node immediately downstream from given node.
"""
function outlets(sn::StreamfallNetwork, node_name::String)::Array{Int}
    nid, _ = sn[node_name]
    return outneighbors(sn.mg, nid)
end


"""
    create_node(mg::MetaDiGraph, node_name::String, details::AbstractDict, nid::Int)

Create a node specified with given name (if it does not exist).

Returns
- `this_id` : ID of node (if pre-existing) and
- `nid` : incremented node id for entire network (equal to `this_id` if exists)
"""
function create_node(mg::MetaDiGraph, node_name::String, details::AbstractDict, nid::Int)
    details = copy(details)

    match = collect(MetaGraphs.filter_vertices(mg, :name, node_name))
    if isempty(match)
        node_type = details["node_type"]

        dtype = eval(Symbol(node_type))
        n = nothing
        try
            n = dtype(node_name, details)
        catch err
            throw(ArgumentError("Unsupported node type: $(node_type) or unknown node details"))
        end

        # Set function for node if specified
        if haskey(details, "func")
            func_spec = details["func"]
            if func_spec isa Function
                func = func_spec
            elseif func_spec isa String
                func = eval(Meta.parse(func_spec))
            end
        else
            func = run_node!
        end

        set_props!(mg, nid, Dict(:name=>node_name,
                                 :node=>n,
                                 :nfunc=>func))

        this_id = nid
    else
        this_id = match[1]
    end

    return this_id
end


"""
    create_network(name::String, network::AbstractDict)::StreamfallNetwork

Create a StreamNetwork from a YAML-derived specification.

# Example
```julia-repl
julia> using OrderedCollections
julia> network_spec = YAML.load_file("example_network.yml"; dicttype=OrderedDict{Any,Any})
julia> sn = create_network("Example Network", network_spec)
```
"""
function create_network(name::String, network::AbstractDict)::StreamfallNetwork
    num_nodes = length(network)
    mg = MetaDiGraph(num_nodes)
    MetaGraphs.set_prop!(mg, :description, name)

    for (nid, (node, details)) in enumerate(network)
        n_name = string(node)

        this_id = create_node(mg, n_name, details, nid)

        if haskey(details, "inlets")
            inlets = details["inlets"]
            if !isnothing(inlets)
                for inlet in inlets

                    inlet_id = findall(keys(network) .== inlet)[1]

                    in_id = create_node(mg, string(inlet), network[inlet], inlet_id)
                    add_edge!(mg, inlet_id => nid)
                end
            end
        end

        if haskey(details, "outlets")
            outlets = details["outlets"]
            if !isnothing(outlets)
                msg = "Streamfall currently only supports a single outlet. ($(length(outlets)))"
                @assert length(outlets) <= 1 || throw(ArgumentError(msg))

                for outlet in outlets
                    outlet_id = findall(keys(network) .== outlet)[1]
                    out_id = create_node(mg, string(outlet), network[outlet], outlet_id)
                    add_edge!(mg, nid => outlet_id)
                end
            end
        end
    end

    sn = StreamfallNetwork(mg)

    return sn
end


"""
    reset!(sn::StreamfallNetwork)::Nothing

Reset a network.
"""
function reset!(sn::StreamfallNetwork)::Nothing
    mg = sn.mg
    v_ids = vertices(mg)
    for i in v_ids
        curr_node = MetaGraphs.get_prop(mg, i, :node)
        reset!(curr_node)
    end
end


function extract_node_spec!(sn::StreamfallNetwork, nid::Int, spec::AbstractDict)::Nothing
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
    network_spec = OrderedDict(
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
function extract_network_spec(sn::StreamfallNetwork)::OrderedDict
    _, outlets = find_inlets_and_outlets(sn)
    spec = OrderedDict()
    for nid in outlets
        extract_node_spec!(sn, nid, spec)
    end

    return spec
end


function save_network_spec(sn::StreamfallNetwork, fn::String)::Nothing
    spec = extract_network_spec(sn)
    write_file(fn, spec)

    return nothing
end


Base.show(io::IO, ::MIME"text/plain", sn::StreamfallNetwork) = show(io, sn)
function Base.show(io::IO, sn::StreamfallNetwork)

    name = MetaGraphs.get_prop(sn.mg, :description)

    println(io, "Network Name: $(name)")
    println(io, "Represented Area: $(area(sn))")
    println("-"^17, "\n")

    vs = vertices(sn.mg)

    show_verts = vs
    if length(vs) > 4
        show_verts = [1, 2, nothing, length(vs)-1, length(vs)]
    end

    for nid in show_verts
        if isnothing(nid)
            println(io, "⋮ \n")
            continue
        end
        println(io, "Node $(nid): \n")
        show(io, sn[nid])
    end
end


"""
    plot_network(sn::StreamfallNetwork)

Simple plot of stream network.
"""
function plot_network(sn::StreamfallNetwork; as_html=false)
    node_labels = ["$(sn[i].name)\n"*string(nameof(typeof(sn[i]))) for i in vertices(sn.mg)]

    if as_html
        plot_func = gplothtml
    else
        plot_func = gplot
    end

    plot_func(sn.mg, nodelabel=node_labels)
end


"""
    save_figure(sn::StreamfallNetwork, fn::String)

Save a figure of the network in SVG format.
"""
function save_figure(sn::StreamfallNetwork, fn::String)
    draw(SVG(fn, 16cm, 16cm), plot_network(sn))
end


function Base.iterate(sn::StreamfallNetwork, state=1)
    if state > length(sn)
        return nothing
    end

    return (sn[state], state+1)
end


Base.length(sn::StreamfallNetwork) = nv(sn.mg)



