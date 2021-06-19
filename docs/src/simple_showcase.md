# A simple example

In this example we showcase a two-node network to represent water levels at Lake Eppalock, 
a dam in the Lower Campaspe catchment located in North-Central Victoria, Australia.

The map below shows Lake Eppalock, along with relevant gauge locations/data 
(click the brown dots on the map to see further gauge details).

This example uses the setup as detailed in [Calibration setup](@ref).

```@raw html
<iframe style="width: 720px; height: 600px; border: none;" src="https://nationalmap.gov.au/#share=s-dIbct7mdo25m7ZK2EVr7Koi4cMp" allowFullScreen mozAllowFullScreen webkitAllowFullScreen></iframe>
```


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
