# A simple showcase of a hydrological system

Here we showcase a two-node network representing inflows into a dam and the dam itself.

The Lower Campaspe catchment - a small semi-arid basin in North-Central Victoria, Australia - is used for the example.

```@raw html
<iframe style="width: 720px; height: 600px; border: none;" src="https://nationalmap.gov.au/#share=s-dIbct7mdo25m7ZK2EVr7Koi4cMp" allowFullScreen mozAllowFullScreen webkitAllowFullScreen></iframe>
```

In this example, we are focused on representing dam levels.

This example uses the setup as detailed in [Calibration setup](@ref)

```julia
@info "Running example stream..."

reset!(sn) # clear any previous runs

# Run the dam node and above
dam_id, dam_node = get_gauge(sn, "406000")
run_node!(sn, dam_id, climate; water_order=hist_dam_releases)

# Get performance metrics
h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

rmse_score = Streamfall.RMSE(h_data, n_data)
nnse_score = Streamfall.NNSE(h_data, n_data)
nse_score = Streamfall.NSE(h_data, n_data)

rmse = round(rmse_score, digits=4)
nnse = round(nnse_score, digits=4)
nse = round(nse_score, digits=4)

@info "Scores:" rmse_score nnse_score nse_score


# Results of model run
plot(h_data,
     legend=:bottomleft,
     title="Calibrated IHACRES\n(RMSE: $(rmse); NSE: $(nse))",
     label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

plot!(n_data, label="IHACRES")
```

![](assets/calibrated_example.png)
