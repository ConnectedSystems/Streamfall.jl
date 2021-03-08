# Streamfall.jl

Experimental graph-based streamflow modelling system written in Julialang.

Currently the only model available is the IHACRES rainfall-runoff model, leveraging [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim).

Stream networks are specified in YAML files, with connectivity defined as a single item or a list of entries:

```yaml
Example_node3:
    inlets:
        - Example_node1
        - Example_node2
    outlets: Example_node4
```

Each node definition can then hold the relevant parameter details/values for a given node.

A full example of the spec is available [here](https://github.com/ConnectedSystems/Streamfall.jl/blob/main/test/data/campaspe/campaspe_network.yml).

[LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used for network traversal/analysis.

In the (near) future, it will be possible to read in stream network information from a standardised XML-based format and/or a Shapefile.


## Running a network

The typical use pattern is to identify the outlets for a given network...

```julia
inlets, outlets = find_inlets_and_outlets(g)
```

... and call `run_node!` for each outlet (with relevant climate data), which will recurse through all relevant nodes.

```julia
@info "Running example stream..."
timesteps = sim_length(climate)
for ts in (1:timesteps)
    for outlet in outlets
        run_node!(mg, g, outlet, climate, ts)
    end
end
```

Preliminary usage examples are provided in the [examples](https://github.com/ConnectedSystems/Streamfall.jl/tree/main/examples) directory.
