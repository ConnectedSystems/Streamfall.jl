<div align="center">
<img align="center" src="docs/src/assets/logo.svg" alt="Streamfall.jl" fill="currentColor"/>
<p>A graph-based streamflow modelling system written in Julialang.</p>

[![DOI](https://zenodo.org/badge/345341654.svg)](https://zenodo.org/badge/latestdoi/345341654)  [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://connectedsystems.github.io/Streamfall.jl/dev)

</div>

Streamfall leverages the Julia language and ecosystem to provide:
- Quick heterogeneous modelling of a stream network
- Use of different rainfall-runoff models and their ensembles in tandem
- Support for modelling and assessment of interacting systems
- A wide range of performance metrics

This package includes implementations of the following:
- GR4J
- HyMod
- IHACRES
- SIMHYD

Performance is expected to be similar to implementations in C and Fortran.

Naive timings (using `@time`) for an example dataset spanning 1963-07-05 - 2014-12-31 (18808 days, approximately 51.5 years)

- SimpleHyModNode \
  0.016502 seconds (469.25 k allocations: 12.902 MiB)
- GR4JNode \
  0.015274 seconds (224.75 k allocations: 5.584 MiB)
- SIMHYDNode \
  0.039540 seconds (638.01 k allocations: 15.190 MiB, 46.99% gc time)
- IHACRESBilinearNode \
  0.021734 seconds (675.63 k allocations: 17.773 MiB)

The IHACRES rainfall-runoff model was previously implemented with [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim) but has since been ported to pure Julia.

[Graphs](https://github.com/JuliaGraphs/Graphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used underneath for network traversal/analysis.

## Installation

Streamfall is now registered! The latest release version can be installed with:

```julia
] add Streamfall
```

or the latest development version from GitHub with `dev` () or `add`:

```julia
# Editable install
] dev Streamfall#main

] add https://github.com/ConnectedSystems/Streamfall.jl#main
```

## Development

Local development should follow the usual process of git cloning the repository.

To build locally:

```bash
$ julia --project=.
julia>] build
```

To run tests:

```bash
julia>] test
```

## Usage

The examples below use data from the CAMEL-AUS dataset, available here:

> Fowler, K. J. A., Acharya, S. C., Addor, N., Chou, C., and Peel, M. C.: CAMELS-AUS: hydrometeorological time   series and landscape attributes for 222 catchments in Australia, Earth Syst. Sci. Data, 13, 3847–3867, https://doi.org/10.5194/essd-13-3847-2021, 2021.

Note that since start of development, an updated dataset is incoming (currently under review):

> Fowler, K. J. A., Zhang, Z., and Hou, X.: CAMELS-AUS v2: updated hydrometeorological timeseries and landscape attributes for an enlarged set of catchments in Australia, Earth Syst. Sci. Data Discuss. [preprint], https://doi.org/10.5194/essd-2024-263, in review, 2024.

Climate data was sourced from the Climate [Change in Australia](https://www.climatechangeinaustralia.gov.au)
data service. Additional data was extracted from the [Long Paddock data silo](https://www.longpaddock.qld.gov.au/silo/).

## Quick start (single node)

The examples below are run from the [examples](https://github.com/ConnectedSystems/Streamfall.jl/tree/main/examples) directory.

```julia
using Statistics
using CSV, DataFrames, YAML
using StatsPlots
using Streamfall

data_dir = joinpath(dirname(dirname(pathof(Streamfall))), "test/data")

# Load data file which holds observed streamflow, precipitation and PET data
obs_data = CSV.read(
    joinpath(data_dir, "cotter/climate/CAMELS-AUS_410730.csv"), 
    DataFrame; 
    comment="#"
)
# 18808×8 DataFrame
#    Row │ year   month  day    Date        410730_P    410730_PET  410730_max_T  410730_Q
#        │ Int64  Int64  Int64  Date        Float64     Float64     Float64       Float64
# ───────┼─────────────────────────────────────────────────────────────────────────────────
#      1 │  1963      7      5  1963-07-05   0.204475     1.02646        6.80409  127.322
#      2 │  1963      7      6  1963-07-06   4.24377      0.790078       5.91556  110.224
#      3 │  1963      7      7  1963-07-07   5.20097      0.400584       3.02218  117.653

# By default, Streamfall expects Date, Precipitation (P), Evapotranspiration (ET) and
# Flow (Q) columns for each gauge as a minimum.
# Some rainfall-runoff models may also require Temperature (T) data.

Qo = extract_flow(obs_data, "410730")
climate = extract_climate(obs_data)

# Create a node (node type, node id, catchment/watershed area)
hymod_node = create_node(SimpleHyModNode, "410730", 129.2);

# Attempt to fit model parameters for 30 seconds
# Here, we use RMSE (Root Mean Square Error) as the objective function
# Note that Streamfall assumes any objective function is to be minimized.
calibrate!(hymod_node, climate, Qo, Streamfall.RMSE; MaxTime=30)

# Run the fitted model
run_node!(hymod_node, climate)

# Display a basic overview plot (shows time series and Q-Q plot)
# using a 366 day offset (e.g., ~1 year burn-in period)
quickplot(Qo, hymod_node, climate, "HyMod"; burn_in=366)

# Save figure
savefig("quick_example.png")
```

## Quick start (network of nodes)

```julia
# Load and generate stream network
network_spec = YAML.load_file("network.yml")
sn = create_network("Example Network", network_spec)

# Show figure of network
plot_network(sn)

# Calibrate network using the BlackBoxOptim package
# keyword arguments will be passed to the `bboptimize()` function
calibrate!(sn, climate, Qo, Streamfall.RMSE; MaxTime=180.0)

# Run stream network
# There is also `run_catchment!()` which does the same thing
run_basin!(sn, climate)

# Get a specific node in network
node = sn[1]  # get the first node in the network ("node1")

# Nodes can also be retrieved by name
# which will also return its position in the network:
# nid, node = sn["node1"]

# Compare "goodness-of-fit"
Streamfall.RMSE(obs_streamflow, node.outflow)

# Save calibrated network spec to file
Streamfall.save_network(sn, "calibrated_example.yml")
```

To display an overview of a node or network:

```julia
julia> node
Name: 406219 [IHACRESBilinearNode]
Area: 1985.73
┌──────────────┬───────────┬─────────────┬─────────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│    Parameter │     Value │ Lower Bound │ Upper Bound │                                                                                                              Description │
├──────────────┼───────────┼─────────────┼─────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│            d │   84.2802 │        10.0 │       550.0 │ Catchment moisture deficit threshold, higher values indicate the catchment can hold more water before generating runoff. │
│           d2 │   2.42241 │      0.0001 │        10.0 │               Scaling factor (d*d2) which creates a second threshold, changing the shape of effective rainfall response. │
│            e │  0.812959 │         0.1 │         1.5 │                      PET conversion factor, controls the rate of evapotranspiration losses, converts temperature to PET. │
│            f │   2.57928 │        0.01 │         3.0 │                             Plant stress threshold, controls at what moisture deficit plants begin to experience stress. │
│            a │   5.92338 │         0.1 │        10.0 │                                    Quickflow storage coefficient, where higher values lead to faster quickflow response. │
│            b │ 0.0989926 │       0.001 │         0.1 │                                            Slowflow storage coefficient, lower values lead to slower baseflow recession. │
│ storage_coef │   1.86134 │     1.0e-10 │        10.0 │                               Groundwater interaction factor, controling how water is exchanged with deeper groundwater. │
│        alpha │  0.727905 │      1.0e-5 │         1.0 │                                                      Effective rainfall scaling factor, partitions rainfall into runoff. │
└──────────────┴───────────┴─────────────┴─────────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### Network specification

Stream networks are specified as Dictionaries, with an entry for each node.

An example spec from a YAML file is shown here, with connectivity between nodes defined by their names.

```yaml

# ... Partial snippet of stream definition as an example ...

Node3:
    node_type: IHACRESBilinearNode  # node type, typically tied to the attached model
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

Each node defines a subcatchment and holds the relevant parameter values for the associated
model.

## Running a network

```julia
using CSV, DataFrames, YAML
using StatsPlots
using Streamfall


# Load a network from a file, providing a name for the network and the file path.
# Creates a graph representation of the stream with associated metadata.
data_dir = joinpath(dirname(dirname(pathof(Streamfall))), "test/data")
sn = load_network("Example Network", joinpath(data_dir, "campaspe/campaspe_network.yml"))

# Load climate data, in this case from a CSV file with data for all nodes.
# Indicate which columns are precipitation and evaporation data based on partial identifiers
climate = Climate(joinpath(data_dir, "campaspe/climate/climate.csv"), "_rain", "_evap")

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
inflow = ...      # inflows for each time step
extractions = ... # extractions from stream for each time step
gw_flux = ...     # forced groundwater interactions for each time step
run_node!(node, climate; inflow=inflow, extraction=extractions, exchange=gw_flux)
```

Another approach is to identify the outlets for a given network...

```julia
inlets, outlets = find_inlets_and_outlets(sn)
```

... and call `run_node!` for each outlet (with relevant climate data), which will recurse
through all relevant nodes upstream.


```julia
@info "Running example stream..."
timesteps = sim_length(climate)

reset!(sn)
prep_state!(sn, timesteps)
for ts in (1:timesteps)
    for outlet in outlets
        run_node!(sn, outlet, climate, ts)
    end
end
```

See the [docs](https://connectedsystems.github.io/Streamfall.jl/dev) for an overview and example applications.

Further preliminary usage examples are provided in the [examples](https://github.com/ConnectedSystems/Streamfall.jl/tree/main/examples) directory.
