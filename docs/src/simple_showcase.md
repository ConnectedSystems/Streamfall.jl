# A simple example

In this example we showcase a two-node network: A gauge providing inflows into a dam.

The dam in question is Lake Eppalock, in the Lower Campaspe catchment located in
North-Central Victoria, Australia.

The map below shows Lake Eppalock, along with relevant gauge locations/data
(click the markers on the map to see further gauge details).

This example uses the results as detailed in [Calibration setup](@ref).

```@raw html
<iframe style="width: 720px; height: 600px; border: none;" src="https://nationalmap.gov.au/#share=s-kxvHElDvlHdB4D4XslDCT70YHZ3" allowFullScreen mozAllowFullScreen webkitAllowFullScreen></iframe>
```

```julia
"""
This script is run in the `examples` directory.
"""

using CSV, DataFrames, YAML
using Plots
using Streamfall

# Load climate data - in this case from a CSV file with data for all nodes.
climate_data = CSV.read(
    "../test/data/campaspe/climate/climate.csv",
    DataFrame;
    comment="#"
)

# Indicate which columns are precipitation and evaporation data based on partial identifiers
climate = Climate(climate_data, "_rain", "_evap")

calib_data = CSV.read(
    "../test/data/campaspe/gauges/outflow_and_level.csv",
    DataFrame;
    comment="#"
)

# Historic extractions from the dam
extraction_data = CSV.read("../test/data/campaspe/gauges/dam_extraction.csv", DataFrame; comment="#")

# Load the two-node example network
sn = load_network("Example Network", "calibration/lake_eppalock.yml")

# Run the dam node and above
dam_id, dam_node = sn["406000"]
run_node!(sn, dam_id, climate; extraction=extraction_data)

# Get performance metrics
dam_obs = calib_data[:, "406000"]
dam_sim = dam_node.level

rmse_score = Streamfall.RMSE(dam_obs, dam_sim)
nnse_score = Streamfall.NNSE(dam_obs, dam_sim)
nse_score = Streamfall.NSE(dam_obs, dam_sim)

rmse = round(rmse_score, digits=4)
nnse = round(nnse_score, digits=4)
nse = round(nse_score, digits=4)

@info "Scores:" rmse nnse nse

# Results of model run
quickplot(dam_obs, dam_sim, climate, "IHACRES", false; burn_in=366)
```

![](assets/calibrated_example.png)
