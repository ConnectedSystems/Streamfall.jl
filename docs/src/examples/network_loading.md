# Network Loading

Loading a pre-defined network from a YAML file.

```julia
using YAML
using StatsPlots, GraphPlot
using Streamfall

network_file = joinpath(
    dirname(dirname(pathof(Streamfall))),
    "test/data/campaspe/campaspe_network.yml"
)

sn = load_network("Example Network", network_file)

# Find all inlets and outlets
inlets, outlets = find_inlets_and_outlets(sn)

# Show figure of network
plot_network(sn)
```