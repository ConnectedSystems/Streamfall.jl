# Node creation

Example of creating a node with default parameters, then updating the parameters with
user-defined values.

```julia
using Streamfall

hymod_node = create_node(SimpleHyModNode, "410730", 129.2)

# Hymod parameters ("hy_" prefix is simply to avoid any variable name conflicts)
hy_Sm_max = 250.0
hy_B = 1.0
hy_alpha = 0.2
hy_Kf = 0.5
hy_Ks = 0.05

# Update parameters
update_params!(hymod_node, hy_Sm_max, hy_B, hy_alpha, hy_Kf, hy_Ks)
```
