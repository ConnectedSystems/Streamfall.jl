<img align="center" src="docs/src/assets/logo.png" alt="Streamfall.jl" />

Streamfall: An experimental graph-based streamflow modelling system written in Julialang.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://connectedsystems.github.io/Streamfall.jl/dev)

Aims of the project are to leverage the Julia language and ecosystem to support:
- Quick application and exploratory analysis
- Use of different rainfall-runoff models in tandem [**aspiration**]
- Modelling and assessment of interacting systems
- Parallel scenario runs

**Note:** the only model currently available is the IHACRES rainfall-runoff model, leveraging [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim).

[LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used underneath for network traversal/analysis.


Development version of the documentation can be found [here](https://connectedsystems.github.io/Streamfall.jl/dev).


Stream networks are specified in YAML files, with connectivity defined as a single item or a list of entries:

```yaml

# ... Partial snippet of stream definition as an example ...

Node3:
    node_type: IHACRESNode  # node type, typically tied to the attached model
    inlets:  # nodes that contribute incoming streamflow
        - Node1
        - Node2
    outlets: Node4  # node that this node flows to
    area: 150.0  # subcatchment area in km^2
    parameters:
        # model specific parameters defined here
        ...
```

A full example of the spec is available [here](https://github.com/ConnectedSystems/Streamfall.jl/blob/main/test/data/campaspe/campaspe_network.yml). The snippet above defines `Node 3` in the diagram below.

[![](https://mermaid.ink/img/eyJjb2RlIjoiZ3JhcGggTFJcbiAgICBBKChOb2RlIDEpKSAtLT4gQygoTm9kZSAzKSlcbiAgICBCKChOb2RlIDIpKSAtLT4gQ1xuICAgIEMgLS0-IEQoKE5vZGUgNCkpXG4gICAgXG4gICIsIm1lcm1haWQiOnsidGhlbWUiOiJkZWZhdWx0In0sInVwZGF0ZUVkaXRvciI6ZmFsc2UsImF1dG9TeW5jIjp0cnVlLCJ1cGRhdGVEaWFncmFtIjpmYWxzZX0)](https://mermaid-js.github.io/mermaid-live-editor/edit##eyJjb2RlIjoiZ3JhcGggTFJcbiAgICBBKChOb2RlIDEpKSAtLT4gQygoTm9kZSAzKSlcbiAgICBCKChOb2RlIDIpKSAtLT4gQ1xuICAgIEMgLS0-IEQoKE5vZGUgKSlcbiAgICBcbiAgIiwibWVybWFpZCI6IntcbiAgXCJ0aGVtZVwiOiBcImRlZmF1bHRcIlxufSIsInVwZGF0ZUVkaXRvciI6ZmFsc2UsImF1dG9TeW5jIjp0cnVlLCJ1cGRhdGVEaWFncmFtIjpmYWxzZX0)



Each node defines a subcatchment and holds the relevant parameter values for the associated model. In the future, it will be possible to read in stream network information from other formats (e.g., GeoPackage).


## Running a network



```julia
using YAML, DataFrames, CSV
using Streamfall


# Load and generate stream network
# Creates a graph representation of the stream with associated metadata.
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")

# Name of network/catchment and its specification
sn = create_network("Example Network", network)

# Load climate data - in this case from a CSV file with data for all nodes.
climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv",
                          comment="#",
                          dateformat="YYYY-mm-dd"))

# Indicate which columns are precipitation and evaporation data based on partial identifiers
climate = Climate(climate_data, "_rain", "_evap")


# This runs an entire stream network
@info "Running an example stream..."
run_catchment!(sn, climate)

@info "Displaying outflow from node 406219"
node_id, node = get_gauge(sn, "406219")
plot(node.outflow)
```

Individual nodes can be run for more fine-grain control.

```julia
# Run up to a point in the stream for all time steps.
# All nodes upstream will be run as well (but not those downstream)
node_id, node = get_gauge(sn, "406219")
run_node!(sn, node_id, climate)

# Reset a node (clears stored states)
reset!(node)

# Run a specific node, and only a specific node, for all time steps
inflow = ...      # array of inflows for each time step
extractions = ... # extractions from stream for each time step
gw_flux = ...     # forced groundwater interactions for each time step
run_node!(node, climate; inflow=inflow, water_order=extractions, exchange=gw_flux)
```

Another approach is to identify the outlets for a given network...

```julia
inlets, outlets = find_inlets_and_outlets(sn)
```

... and call `run_node!` for each outlet (with relevant climate data), which will recurse through all relevant nodes upstream.


```julia
@info "Running example stream..."
timesteps = sim_length(climate)
for ts in (1:timesteps)
    for outlet in outlets
        run_node!(sn, outlet, climate, ts)
    end
end
```

Preliminary usage examples are provided in the [examples](https://github.com/ConnectedSystems/Streamfall.jl/tree/main/examples) directory.
