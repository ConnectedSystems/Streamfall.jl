# Node creation

Example of creating a node with default parameters, then updating the parameters with
user-defined values.

```julia
using Streamfall

hymod_node = create_node(SimpleHyModNode, "410730", 129.2)

# Hymod parameters
Sm_max = 250.0
B = 1.0
alpha = 0.2
Kf = 0.5
Ks = 0.05

# Update parameters
update_params!(hymod_node, Sm_max, B, alpha, Kf, Ks)
```
