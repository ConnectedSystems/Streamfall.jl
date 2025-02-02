# Revisiting weighted ensembles

In the previous section, a simple approach to developing ensembles of models is introduced.

Findings in Arsenault et al. (2015) and Wan et al. (2021) suggest that combinations of
models and objective functions are better than a single model and objective, even with a
simple arithmetic averaging approach.

While in the previous section two model types are applied, they both were calibrated with
the same objective function. Here, we revisit the example applying the non-parametric
formulation of the Kling-Gupta Efficiency metric (npKGE) as it is regarded as being better
suited for low-flow conditions to one of the models in the ensemble.

An alternative could involve Box-Cox transformed flows.

```julia
using CSV, Plots
using DataFrames
using Streamfall

climate = Climate("../test/data/campaspe/climate/climate.csv", "_rain", "_evap")

# Historic flows and dam level data
obs_data = CSV.read(
    "../test/data/cotter/climate/CAMELS-AUS_410730.csv",
    DataFrame;
    comment="#"
)

Qo = extract_flow(obs_data, "410730")
climate = extract_climate(obs_data)

# Create one instance each of IHACRES_CMD and GR4J
ihacres_node = create_node(IHACRESBilinearNode, "410730", 129.2)
gr4j_node = create_node(GR4JNode, "410730", 129.2)
symhyd_node = create_node(SYMHYDNode, "410730", 129.2)
hymod_node = create_node(SimpleHyModNode, "410730", 129.2)

# Calibrate with different objective functions
# NmKGE for IHACRES and NnpKGE for GR4J
calibrate!(ihacres_node, climate, Qo, (obs, sim) -> 1.0 .- Streamfall.NmKGE(obs, sim); MaxTime=180)
calibrate!(gr4j_node, climate, Qo, (obs, sim) -> 1.0 .- Streamfall.NnpKGE(obs, sim); MaxTime=180)
calibrate!(symhyd_node, climate, Qo, Streamfall.RMSE; MaxTime=180)
calibrate!(hymod_node, climate, Qo, (obs, sim) -> Streamfall.inverse_metric(obs, sim, (obs, sim) -> 1.0 .- Streamfall.NmKGE(obs, sim); comb_method=mean))

# Create an ensemble
ensemble = create_node(GREnsembleNode, [ihacres_node, gr4j_node, symhyd_node, hymod_node])

# Calibrate the ensemble weights
# calibrate!(ensemble, climate, Qo, (obs, sim) -> 1.0 .- Streamfall.NmKGE(obs, sim); MaxTime=180)
calibrate!(ensemble, climate, Qo, Streamfall.RMSE)
```


```julia
# run_node!(ihacres_node, climate)
# run_node!(gr4j_node, climate)
run_node!(ensemble, climate)

burn_in = 365
burn_dates = timesteps(climate)[burn_in:end]
burn_obs = Qo[burn_in:end, "410730"]

ensemble_qp = quickplot(burn_obs, ensemble.outflow[burn_in:end], climate, "GRC Ensemble", true)

ensemble_xs = temporal_cross_section(burn_dates, burn_obs, ensemble.outflow[burn_in:end]; title="GRC Ensemble (IHACRES-GR4J)")

plot(ensemble_qp, ensemble_xs; layout=(2, 1), size=(800, 600))
```

## Additional remarks

The approach to weighting and averaging model outputs also play a role with both identifying
Granger–Ramanathan average variant C (GRC) to be the most performant with respect to the
Nash-Sutcliffe Efficiency metric.

While GRC is not offered in Streamfall it could hypothetically be implemented by the user.

## References

1. Arsenault, R., Gatien, P., Renaud, B., Brissette, F., Martel, J.-L., 2015. \
   A comparative analysis of 9 multi-model averaging approaches in hydrological continuous \
   streamflow simulation. Journal of Hydrology 529, 754–767. \
   https://doi.org/10.1016/j.jhydrol.2015.09.001

2. Wan, Y., Chen, J., Xu, C.-Y., Xie, P., Qi, W., Li, D., Zhang, S., 2021. Performance \
   dependence of multi-model combination methods on hydrological model calibration strategy \
   and ensemble size. Journal of Hydrology 603, 127065. \
   https://doi.org/10.1016/j.jhydrol.2021.127065
