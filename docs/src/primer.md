# Primer

Streamfall is a stream network modelling framework with integrated systems analysis and modelling in mind. The aim is to simplify the construction of basin-scale hydrology models, itself a constituent of a larger system of systems. The overarching concepts are explained here.

Primary components of the Streamfall framework include:

1. The graph network defining the stream and the model associated with each node in the network
2. Data for the basin, including climate data and hydrologic interactions driven by other systems
3. The functions which run the network as a whole and individual nodes

## Defining a network

A stream network is defined through a YAML specification file (the "spec").

A single node network is shown below using a `IHACRES` model (see [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim) for details on the model).

The spec takes the following form:

```YAML
# Node name (the Gauge ID is used here)
410730:
    # The node type which defines which model is used for this node
    # In this case, it is the IHACRES with the bilinear formulation of the CMD module
    node_type: BilinearNode

    # This spec defines a single node system
    # so it has no nodes upstream (inlets) or downstream (outlets)
    inlets:
    outlets:
    area: 130.0  # subcatchment area in km^2 (from BoM)

    # Model parameters (in this case, for IHACRES)
    parameters:    
        d: 200.0     # millimeters
        d2: 2.0      # multiplier applied to `d`
        e: 1.0       # ET scaling factor, dimensionless
        f: 0.8       # multiplier applied to `d` to determine effective rainfall, dimensionless
        a: 0.9       # quickflow scaling factor
        b: 0.1       # slowflow scaling factor
        storage_coef: 2.9  # groundwater interaction factor
        alpha: 0.95  # effective rainfall scaling factor
        initial_storage: 0.0  # initial CMD value, CMD > 0 means there is a deficit

    # additional node-specific parameters
    # (unused in this example so can be ignored)
    level_params:  
        - -3.3502  # p1
        - 0.68340  # p2
        - 4.50     # p3
        - 5.0      # p4
        - 0.35     # p5
        - 1.41     # p6
        - -1.45    # p7
        - 6.75     # p8
        - 167.845  # CTF
```

The spec is then loaded in Julia and passed into `create_network()`

```julia
# Load file
network = YAML.load_file("network.yml")

# Create network from spec, with a human-readable name.
sn = create_network("Gingera Catchment", network)
```

Printing the network displays a summary of the nodes:

```julia-repl
julia> sn

Network Name: Gingera Catchment
Represented Area: 130.0
-----------------

Node 1 :

Name: 410730 [BilinearNode]
Area: 130.0
┌──────────────┬───────┬─────────────┬─────────────┐
│    Parameter │ Value │ Lower Bound │ Upper Bound │
├──────────────┼───────┼─────────────┼─────────────┤
│            d │ 200.0 │        10.0 │       550.0 │
│           d2 │   2.0 │      0.0001 │        10.0 │
│            e │   1.0 │         0.1 │         1.5 │
│            f │   0.8 │        0.01 │         3.0 │
│            a │   0.9 │         0.1 │        10.0 │
│            b │   0.1 │       0.001 │         0.1 │
│ storage_coef │   2.9 │     1.0e-10 │        10.0 │
│        alpha │  0.95 │      1.0e-5 │         1.0 │
└──────────────┴───────┴─────────────┴─────────────┘
```


Each node will be assigned an internal node identifier based on their order and position
in the network.

Individual nodes can also be created programmatically:

```julia
# Programmatically create a node (from a spec)
new_node = BilinearNode("410730", network["410730"])

# Creating the same node manually by specifying model parameters
# Argument order: node_name, area, d, d2, e, f, a, b, storage_coef, alpha, initial cmd, initial quickflow, initial slowflow, initial gw_store
new_node = BilinearNode("410730", 130.0, 95.578, 1.743, 1.047, 1.315, 99.134, 0.259, 2.9, 0.785, 100.0, 0.0, 0.0, 0.0)
```

Of course, model parameters may not be known in advance.

Streamfall nodes hold parameter information including their usual bounds/ranges.

These can be examined like so:

```julia
param_names, x0, bounds = param_info(node)
```

where `x0` are the current values.

Example calibration approaches are detailed in [Example calibration](@ref)

## Required data

All data is expected in DataFrames with the following convention:

- All time series must have a "Date" column in `YYYY-mm-dd` format.
- By convention, column names follow the pattern: `[node_name]_[phenomenon]_[other metadata]`
- Unit of measure itself is optionally included in square brackets (`[]`)
- In-file comments are indicated with a hash (`#`)


Data that may be optionally provided include:

- `inflow` : specifying the incoming volume of water (when running the model for a specific node)
- `extractions` : additional extractions from the stream
- `exchange` : additional forcing for groundwater interactions.


These may be provided as a Dictionary of arrays with the node name as the key.

More details may be found in the [Input data format](@ref) page.

### Climate data

Climate data is treated as a special case of the above and is required for a scenario to run.
Currently the expectation is that two phenomena are provided for each node (one of which is rainfall).

In this example case, precipitation and evaporation are provided (marked by identifiers `_rain` and `_evap` respectively). Below is an example for a two-node system.

```csv
Date, 406214_rain, 406214_evap, 406219_rain, 406219_evap
1981-01-01, 0.0, 4.8, 0.0, 4.9
1981-01-02, 0.1, 0.5, 0.1, 3.3
1981-01-03, 10.5, 5.3, 7.2, 2.3
1981-01-04, 9.89, 7.9, 6.1, 4.3
1981-01-05, 0.3, 4.2, 0.2, 6.4
... snip ...
```

```julia
# Load data from CSV
date_format = "YYYY-mm-dd"
climate_data = DataFrame!(CSV.File(joinpath(data_path, "climate/climate_historic.csv"),
                          comment="#",
                          dateformat=date_format))

# Create a climate object, specifying which identifiers to use.
climate = Climate(climate_data, "_rain", "_evap")
```

## Running a network or node

To run an entire basin network, without any dynamic interaction with "external" systems:

```julia
run_basin!(sn::StreamfallNetwork, climate::Climate)
```

This will identify the final "outlet" of the stream network and recurse upstream to run all nodes.

Individual nodes can also be run:

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

When interactions with other socio-environmental systems are expected, it can become necessary to run each node individually as needed.

Interactions with these external systems are represented as influencing:

1. Inflows to a node
2. Extraction of water from a stream
3. Additional flux to/from the groundwater system

The following pattern can be used in such a context:


```julia
# ... as an example, this is the 100th day in simulation...
this_timestep = 100

# Inflow comes from upstream
# This could be obtained by running all nodes upstream, for example
#     node_id, node = sn["406219"]
#     run_node!(sn, node_id, climate)
inflow = run_node!(...)

# run external model that provides extraction, say 10ML
# This does not have to be a model in Julia per se.
extractions = some_other_model(...)

# Run a different model which provides groundwater interactions
exchange = another_model(...)

# Obtain outflow and level for this time step
outflow, level = run_node!(target_node, climate, this_timestep; inflow=inflow, extraction=extractions, exchange=exchange)
```

Specific examples can be found in the Examples section.
