# Network Loading

Loading a pre-defined network from a YAML file.

```julia
using YAML
using Streamfall


sn = load_network("Example Network", "../test/data/campaspe/campaspe_network.yml")

# Find all inlets and outlets
inlets, outlets = find_inlets_and_outlets(sn)
```