# Node creation

Example of creating a node with default parameters, then updating the parameters with
user-defined values.

```julia
using Streamfall

hymod_node = create_node(SimpleHyModNode, "410730", 129.2)
# Name: 410730 [SimpleHyModNode]
# Area: 129.2
# ┌───────────┬───────┬─────────────┬─────────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────┐
# │ Parameter │ Value │ Lower Bound │ Upper Bound │                                                                                          Description │
# ├───────────┼───────┼─────────────┼─────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────┤
# │    Sm_max │ 250.0 │         1.0 │       500.0 │                                                                       Maximum soil storage capacity. │
# │         B │   1.0 │         0.0 │         2.0 │                        Controls how quickly the catchment becomes saturated as rainfall accumulates. │
# │     alpha │   0.2 │         0.0 │         1.0 │ The split between quick and slow flow components. Higher values direct more water through quickflow. │
# │        Kf │   0.5 │         0.1 │      0.9999 │                                                                                 Quickflow recession. │
# │        Ks │  0.05 │       0.001 │         0.1 │                                                                                  Slowflow recession. │
# └───────────┴───────┴─────────────┴─────────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────┘

# Hymod parameters ("hy_" prefix is simply to avoid any variable name conflicts)
hy_Sm_max = 370.0
hy_B = 0.5
hy_alpha = 0.3
hy_Kf = 0.25
hy_Ks = 0.25

# Update parameters
update_params!(hymod_node, hy_Sm_max, hy_B, hy_alpha, hy_Kf, hy_Ks)

# The "Value" column indicates model parameters have been updated.
print(hymod_node)
# Name: 410730 [SimpleHyModNode]
# Area: 129.2
# ┌───────────┬───────┬─────────────┬─────────────┬──────────────────────────────────────────────────────────────────────────────────────────────────────┐
# │ Parameter │ Value │ Lower Bound │ Upper Bound │                                                                                          Description │
# ├───────────┼───────┼─────────────┼─────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────┤
# │    Sm_max │ 370.0 │         1.0 │       500.0 │                                                                       Maximum soil storage capacity. │
# │         B │   0.5 │         0.0 │         2.0 │                        Controls how quickly the catchment becomes saturated as rainfall accumulates. │
# │     alpha │   0.3 │         0.0 │         1.0 │ The split between quick and slow flow components. Higher values direct more water through quickflow. │
# │        Kf │  0.25 │         0.1 │      0.9999 │                                                                                 Quickflow recession. │
# │        Ks │  0.25 │       0.001 │         0.1 │                                                                                  Slowflow recession. │
# └───────────┴───────┴─────────────┴─────────────┴──────────────────────────────────────────────────────────────────────────────────────────────────────┘
```
