<img align="center" src="docs/src/assets/logo.png" alt="Streamfall.jl" />

Streamfall: An experimental graph-based streamflow modelling system written in Julialang.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://connectedsystems.github.io/Streamfall.jl/dev)
[![DOI](https://zenodo.org/badge/345341654.svg)](https://zenodo.org/badge/latestdoi/345341654)

Aims of the project are to leverage the Julia language and ecosystem to support:
- Quick application and exploratory analysis
- Use of different rainfall-runoff models and their ensembles in tandem
- Modelling and assessment of interacting systems
- Parallel scenario runs

Streamfall now includes tentative implementations of GR4J, HyMod, IHACRES, and SYMHYD.
The IHACRES rainfall-runoff model is implemented with [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim).

[Graphs](https://github.com/JuliaGraphs/Graphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used underneath for network traversal/analysis.

Development version of the documentation can be found [here](https://connectedsystems.github.io/Streamfall.jl/dev).

> [NOTE] Streamfall is currently in its early stages and under active development. Although it is fairly usable for small networks and single node analyses, things may change drastically and unexpectedly.

## Development

To build locally:

```bash
$ julia --project=.
julia>] build
```

To run tests:

```bash
julia>] test
```


## Quick start (single node)

```julia
using YAML, DataFrames, CSV, Plots
using Statistics
using Streamfall, BlackBoxOptim


using CSV, DataFrames
using Streamfall

# Load data file which holds observed streamflow, precipitation and PET data
obs_data = CSV.read("../test/data/cotter/climate/CAMELS-AUS_410730.csv", DataFrame; comment="#")
# 18808×8 DataFrame
#    Row │ year   month  day    Date        410730_P    410730_PET  410730_max_T  410730_Q
#        │ Int64  Int64  Int64  Date        Float64     Float64     Float64       Float64
# ───────┼─────────────────────────────────────────────────────────────────────────────────
#      1 │  1963      7      5  1963-07-05   0.204475     1.02646        6.80409  127.322
#      2 │  1963      7      6  1963-07-06   4.24377      0.790078       5.91556  110.224
#      3 │  1963      7      7  1963-07-07   5.20097      0.400584       3.02218  117.653

# Streamfall expects Date, Precipitation (P), and flow (Q) columns for each gauge at a minimum.
# Some rainfall-runoff models may also require Temmperature (T) data.

Qo = extract_flow(obs_data, "410730")
climate = extract_climate(obs_data)

# Create a node (node type, node id, catchment/watershed area)
hymod_node = create_node(SimpleHyModNode, "410730", 129.2);

# Attempt to fit model parameters for 30 seconds
calibrate!(hymod_node, climate, Qo; MaxTime=30)

# Run the fitted model
run_node!(hymod_node, climate)

# Display a basic overview plot (shows time series and Q-Q plot)
# using a 365 day offset (e.g., 1 year burn-in period)
quickplot(Qo, hymod_node, climate, "HyMod"; burn_in=366, limit=nothing)

# Save figure
savefig("quick_example.png")
```


## Quick start (network of nodes)

!!! We are aware of bugs in the current implementation. These are being addressed.

```julia
# Load and generate stream network
network_spec = YAML.load_file("network.yml")
sn = create_network("Example Network", network_spec)

# Show figure of network
plot_network(sn)

# Calibrate network using the BlackBoxOptim package
# keyword arguments will be passed to the `bboptimize()` function
calibrate!(sn, climate, Qo; MaxTime=180.0)

# Run stream network
# There is also `run_catchment!()` which does the same thing
run_basin!(sn, climate)

# Get a specific node in network
node = sn[1]  # get node 1

# Could also get node by name
# which will also return its position in the network:
# nid, node = sn["node1"]

# Compare "goodness-of-fit"
Streamfall.RMSE(obs_streamflow, node.outflow)

# Save calibrated network spec to file
Streamfall.save_network_spec(sn, "calibrated_example.yml")
```

To display an overview of a node or network:

```julia
julia> node
Name: 406219 [BilinearNode]
Area: 1985.73
┌──────────────┬───────────┬─────────────┬─────────────┐
│    Parameter │     Value │ Lower Bound │ Upper Bound │
├──────────────┼───────────┼─────────────┼─────────────┤
│            d │   84.2802 │        10.0 │       550.0 │
│           d2 │   2.42241 │      0.0001 │        10.0 │
│            e │  0.812959 │         0.1 │         1.5 │
│            f │   2.57928 │        0.01 │         3.0 │
│            a │   5.92338 │         0.1 │        10.0 │
│            b │ 0.0989926 │       0.001 │         0.1 │
│ storage_coef │   1.86134 │     1.0e-10 │        10.0 │
│        alpha │  0.727905 │      1.0e-5 │         1.0 │
└──────────────┴───────────┴─────────────┴─────────────┘
```

### More information

Stream networks are specified as Dictionaries, with an entry for each node.

An example spec from a YAML file is shown here, with connectivity between nodes defined by their names.

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
using YAML, DataFrames, CSV, Plots
using Streamfall


# Load and generate stream network
# Creates a graph representation of the stream with associated metadata.
network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")

# Name of network/catchment and its specification
sn = create_network("Example Network", network)

# Load climate data - in this case from a CSV file with data for all nodes.
climate_data = CSV.File("../test/data/campaspe/climate/climate_historic.csv",
                        comment="#",
                        dateformat="YYYY-mm-dd") |> DataFrame

# Indicate which columns are precipitation and evaporation data based on partial identifiers
climate = Climate(climate_data, "_rain", "_evap")


# This runs an entire stream network
@info "Running an example stream..."
run_catchment!(sn, climate)

@info "Displaying outflow from node 406219"
node_id, node = sn["406219"]
plot(node.outflow)
```

Individual nodes can be run for more fine-grain control.

```julia
# Run up to a point in the stream for all time steps.
# All nodes upstream will be run as well (but not those downstream)
node_id, node = sn["406219"]
run_node!(sn, node_id, climate)

# Reset a node (clears stored states)
reset!(node)

# Run a specific node, and only a specific node, for all time steps
inflow = ...      # array of inflows for each time step
extractions = ... # extractions from stream for each time step
gw_flux = ...     # forced groundwater interactions for each time step
run_node!(node, climate; inflow=inflow, extraction=extractions, exchange=gw_flux)
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

See the [docs](https://connectedsystems.github.io/Streamfall.jl/dev) for an overview and example applications.


Further preliminary usage examples are provided in the [examples](https://github.com/ConnectedSystems/Streamfall.jl/tree/main/examples) directory.


