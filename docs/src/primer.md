# Primer

Streamfall is a stream network modelling framework with integrated systems analysis and
modelling in mind. The aim is to simplify the construction of basin-scale hydrology models,
itself a constituent of a larger system of systems. The overarching concepts are explained
here.

The motivation for Streamfall is to enable flexible modelling of a hydrological system.\
This includes:

- Representation of catchments with heterogenous combinations of rainfall-runoff models \
  (each sub-catchment may be represented by a different hydrological model)
- Support interaction with other models which may represent groundwater interactions or \
  anthropogenic activity (e.g., water extractions)
- High performance relative to available implementations in R and Python

Primary components of the Streamfall framework include:

1. The graph representing a network of gauges and the associated model
2. Data for the basin, including climate data and hydrologic interactions \
   driven by other systems
3. The functions which run the network as a whole and individual nodes

## Defining a network

A stream network is defined through a YAML specification file (the "spec").

A single node network is shown below using a `IHACRES` model.

The spec takes the following form:

```YAML
# Node name (the Gauge ID is used here)
410730:
    # The node type which defines which model is used for this node
    # In this case, it is the IHACRES with the bilinear formulation of the CMD module
    node_type: IHACRESBilinearNode
    area: 130.0  # subcatchment area in km^2 (from the Australian Bureau of Meteorology)

    # This spec defines a single node system
    # so it has no nodes upstream (inlets) or downstream (outlets)
    inlets:
    outlets:

    # Initial CMD state, CMD > 0 means there is a deficit
    initial_storage: 0.0

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
```

The spec is then loaded with `load_network()`

```julia
# Load network from a spec file, with a human-readable name.
sn = load_network("Gingera Catchment", "network.yml")
```

Printing the network displays a summary of the nodes:

```julia-repl
julia> sn

Network Name: Gingera Catchment
Represented Area: 130.0

Node 1
--------
Name: 410730 [IHACRESBilinearNode]
Area: 130.0
┌──────────────┬───────┬─────────────┬─────────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│    Parameter │ Value │ Lower Bound │ Upper Bound │                                                                                                              Description │
├──────────────┼───────┼─────────────┼─────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│            d │ 200.0 │        10.0 │       550.0 │ Catchment moisture deficit threshold, higher values indicate the catchment can hold more water before generating runoff. │
│           d2 │   2.0 │      0.0001 │        10.0 │               Scaling factor (d*d2) which creates a second threshold, changing the shape of effective rainfall response. │
│            e │   1.0 │         0.1 │         1.5 │                      PET conversion factor, controls the rate of evapotranspiration losses, converts temperature to PET. │
│            f │   0.8 │        0.01 │         3.0 │                             Plant stress threshold, controls at what moisture deficit plants begin to experience stress. │
│            a │   0.9 │         0.1 │        10.0 │                                    Quickflow storage coefficient, where higher values lead to faster quickflow response. │
│            b │   0.1 │       0.001 │         0.1 │                                            Slowflow storage coefficient, lower values lead to slower baseflow recession. │
│ storage_coef │   2.9 │     1.0e-10 │        10.0 │                              Groundwater interaction factor, controlling how water is exchanged with deeper groundwater. │
│        alpha │  0.95 │      1.0e-5 │         1.0 │                                                      Effective rainfall scaling factor, partitions rainfall into runoff. │
└──────────────┴───────┴─────────────┴─────────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```


Each node will be assigned an internal node identifier based on their order and position
in the network.

Individual nodes can also be created programmatically:

```julia
# Programmatically create a node (from a spec)
new_node = IHACRESBilinearNode("410730", network["410730"])

# Creating the same node manually by specifying model parameters
# Argument order: node_name, area, d, d2, e, f, a, b, storage_coef, alpha, initial cmd, initial quickflow, initial slowflow, initial gw_store
new_node = IHACRESBilinearNode("410730", 130.0, 95.578, 1.743, 1.047, 1.315, 99.134, 0.259, 2.9, 0.785, 100.0, 0.0, 0.0, 0.0)
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
# Create a climate object, specifying which identifiers to use.
example_data_dir = joinpath(dirname(dirname(pathof(Streamfall))), "test/data")
climate = Climate(
    joinpath(example_data_dir, "campaspe/climate/climate.csv"), 
    "_rain", "_evap"
)
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

Another, more direct approach, is to identify all outlets for a given network and call
`run_node!()` for each outlet with relevant climate data for each timestep.
All relevant upstream nodes will also be run.

```julia
inlets, outlets = find_inlets_and_outlets(sn)

@info "Running example stream..."
timesteps = sim_length(climate)
prep_state!(sn, timesteps)
for ts in (1:timesteps)
    for outlet in outlets
        run_node!(sn, outlet, climate, ts)
    end
end
```

When interactions with other socio-environmental systems are expected, it can become
necessary to run each node individually as needed.

Interactions with these external systems are represented as influencing:

1. Inflows to a node
2. Extraction of water from a stream
3. Additional flux to/from the groundwater system

The following pattern can be used in such a context:


```julia
@info "Running example stream..."
steps = sim_length(climate)
prep_state!(sn, steps)
for ts in 1:steps
    # Run external model that provides extraction **in the same units**
    # This does not have to be a model in Julia, but inter-language interoperability is
    # outside the scope of this example.
    extractions = some_water_extraction_model(...)

    # Run a different model which provides groundwater interactions
    exchange = a_groundwater_model(...)

    for outlet in outlets
        run_node!(sn, outlet, climate, ts; extraction=extractions, exchange=exchange)
    end
end
```

Specific examples can be found in the Examples section.
