using Streamfall
using Statistics, CSV, DataFrames
import .Analysis: TemporalCrossSection
using Plots

using BenchmarkTools


# load observations
HERE = @__DIR__
DATA_PATH = joinpath(HERE, "../test/data/cotter/")

date_format = "YYYY-mm-dd"
obs_data = CSV.File(joinpath(DATA_PATH, "climate/CAMELS-AUS_410730.csv"),
                    comment="#",
                    dateformat=date_format) |> DataFrame

hist_streamflow = obs_data[:, ["Date", "410730_Q"]]
hist_streamflow = rename(hist_streamflow, "410730_Q"=>"410730")

climate_data = obs_data[:, ["Date", "410730_P", "410730_PET"]]
climate = Climate(climate_data, "_P", "_PET")

obs = hist_streamflow[:, "410730"]
xsect_res = Streamfall.Analysis.TemporalCrossSection(climate_data.Date, obs)
xsect = Streamfall.Analysis.cross_section(xsect_res)


ihacres_node = create_node(BilinearNode, "410730", 100.0)
burn_in = 365

@benchmark run_node!(ihacres_node, climate)

# metric = (obs, sim) -> 1.0 - Streamfall.NmKGE(log.(obs[burn_in:end]), log.(sim[burn_in:end]))

# res, opt = calibrate!(ihacres_node, climate, hist_streamflow; metric=metric, MaxTime=900.0)
# run_node!(ihacres_node, climate)

# h_dates = climate_data.Date
# burn_dates = h_dates[burn_in:end]
# burn_obs = obs[burn_in:end]
# fig = Streamfall.temporal_cross_section(burn_dates, burn_obs; yscale=:log10)

# plot(fig, size=(500, 400))


# mod_xsect_res = TemporalCrossSection(burn_dates, burn_obs, ihacres_node.outflow[burn_in:end])
# offset_vals = mod_xsect_res.offsets

# fig2 = Streamfall.temporal_cross_section(burn_dates, burn_obs, ihacres_node.outflow[burn_in:end]; yscale=:log10, title="NmKGE")

# function postprocess(x, ω)
#     idx = (x .+ ω) .>= 0.0
#     x[idx] = x[idx] .+ ω[idx]

#     idx2 = (x .+ ω) .< 0
#     x[idx2] .= 0.0
#     return x
# end

# fig3 = Streamfall.temporal_cross_section(burn_dates, burn_obs, postprocess(ihacres_node.outflow[burn_in:end], offset_vals); yscale=:log10, title="Offset (>= 0.0)")

# plot(fig2, fig3, size=(1200, 500), link=:y, layout=(1,2), dpi=300)
# savefig("xsection_overview_log.png")

# off_vals = postprocess(ihacres_node.outflow[burn_in:end], offset_vals)
# orig_score = round(Streamfall.NmKGE(log.(burn_obs), log.(ihacres_node.outflow[burn_in:end])), digits=2)
# off_score = round(Streamfall.NmKGE(log.(burn_obs), log.(off_vals)), digits=2)

# plot(burn_dates, Streamfall.ME.(burn_obs, ihacres_node.outflow[burn_in:end]), alpha=0.4, label="IHACRES (log NmKGE: $(orig_score))")
# plot!(burn_dates, Streamfall.ME.(burn_obs, off_vals), alpha=0.4, label="Post-Processed (log NmKGE: $(off_score))", dpi=300)
# savefig("t_offset_streamflow_ME_log.png")

