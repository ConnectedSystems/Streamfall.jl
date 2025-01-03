using Plots

using Statistics
import Dates: month, monthday, week, dayofyear
using DataFrames, CSV, Serialization

using Streamfall
import Streamfall: TemporalCrossSection, BlackBoxOptim


HERE = @__DIR__
DATA_PATH = joinpath(HERE, "../test/data/cotter/")

# Load observations
date_format = "YYYY-mm-dd"
obs_data = CSV.File(joinpath(DATA_PATH, "climate/CAMELS-AUS_410730.csv"),
    comment="#",
    dateformat=date_format) |> DataFrame

hist_streamflow = obs_data[:, ["Date", "410730_Q"]]
climate_data = obs_data[:, ["Date", "410730_P", "410730_PET"]]
climate = Climate(climate_data, "_P", "_PET")

burn_in = 366  # 1 year burn-in period + 1 day
obs = hist_streamflow[:, "410730_Q"]
log_obs = log.(obs .+ 0.02)

if !isdefined(Main, :bl_overview) || isnothing(bl_overview)
    if !isfile("uncert_baseline_opt.node")
        # Create objective function to minimize (here we use Normalized KGE')
        func = (obs, sim) -> 1.0 - Streamfall.NmKGE(obs[burn_in:end], sim[burn_in:end])

        # Create an individual node
        ihacres_node = create_node(IHACRESBilinearNode, "410730", 129.2)

        # Calibrate node
        res, opt = calibrate!(ihacres_node, climate, obs, func; MaxTime=60)
        update_params!(ihacres_node, best_candidate(res)...)
        reset!(ihacres_node)
        run_node!(ihacres_node, climate)

        bl_sim = ihacres_node.outflow

        # split_NmKGE_results = Streamfall.naive_split_metric(h_data, n_data; n_members=365, metric=Streamfall.NmKGE)
        # bs_NmKGE = (x) -> Streamfall.NmKGE(x[:, 1], x[:, 2])

        open("uncert_baseline_opt.node", "w") do fh
            serialize(fh, ihacres_node)
        end
    else
        ihacres_node = open("uncert_baseline_opt.node", "r") do fh
            deserialize(fh)
        end

        bl_sim = copy(ihacres_node.outflow)
    end
end

bl_xsect = Streamfall.temporal_cross_section(climate_data.Date, obs, bl_sim; show_extremes=false, yscale=:log10)
bl_overview = quickplot(obs, ihacres_node, climate, "NmKGE", true)

# (x_section, min_section, med_section, max_section, high_uncert, cert, std_error, weights) = temporal_uncertainty(climate_data.Date, obs, bl_sim; period=monthday, threshold=0.3)
# (x_section, min_section, max_section, whisker_range, weights, rough) = temporal_uncertainty(climate_data.Date, obs; threshold=0.3)

# function std_scale(x)
#     return (x .- mean(x)) ./ std(x)
# end

# function weighted_obj_function(obs, sim; metric=Streamfall.NmKGE)
#     x_section = Streamfall.TemporalCrossSection(climate_data.Date[burn_in:end], obs, sim)

#     n = length(obs)
#     weights = 1.0
#     weighted_res = [metric([obs[i]], [sim[i]]) for i in 1:n]
#     weighted_res = weighted_res .* weights

#     mean_weighted_res = mean(weighted_res)
#     # cv = std(weighted_res) / mean_weighted_res

#     scaled = mean_weighted_res
#     # scaled = mean((w0 .- minimum(w0)) / (maximum(w0) - minimum(w0)))
#     # return (scaled*0.9 + (1.0 - ℯ^-cv)*0.075 + (1.0 - ℯ^-rough)*0.025)

#     # scaled_avg_x = 1.0 - median(1.0 ./ abs.(x_section))
#     # scaled_x = median(abs.(x_section))
#     # scaled_med_x = median(minmax_scale(abs.(x_section)))
#     # scaled_avg_min = 1 - mean(minmax_scale(min_section))  # (1.0 - ℯ^-abs(mean(min_section)))
#     # scaled_avg_max = mean(minmax_scale(max_section))
#     upp_section = x_section.cross_section[:, :upper_95]
#     low_section = x_section.cross_section[:, :lower_95]

#     inner_range = upp_section .- low_section
#     cv_range = std(inner_range) / mean(inner_range)
#     # max_range = mean(max_section .- min_section)
#     # min_range = -mean(low_section)
#     # max_range = mean(upp_section)

#     score = scaled # + cv_range  # + (1.0 - ℯ^-rough)  # + min_range + max_range

#     return score
# end

plot(bl_xsect)
savefig("baseline_xsection.png")

# inverse_func(obs, sim) = Streamfall.mean_NmKGE(obs, sim)
# split_func(obs, sim) = Streamfall.naive_split_metric(obs, sim; metric=Streamfall.NmKGE, n_members=14)
function v_metric(obs, sim; metric=Streamfall.NmKGE, n_members=180)
    chunked_scores = Streamfall.naive_split_metric(obs, sim, n_members, metric)
    return mean(chunked_scores)*0.9 + std(chunked_scores)*0.1
end

uncert_objfunc = (obs, sim) -> 1.0 - v_metric(obs, sim)

reset!(ihacres_node)
res, opt = calibrate!(ihacres_node, climate, obs, uncert_objfunc; MaxTime=180)

reset!(ihacres_node)
run_node!(ihacres_node, climate)

# log_sim = log.(ihacres_node.outflow .+ 0.02)
uc_xsect = Streamfall.temporal_cross_section(climate_data.Date, obs, ihacres_node.outflow; show_extremes=false, yscale=:log10)
plot(bl_xsect, uc_xsect, link=:y, layout=(1, 2), size=(850, 400))
savefig("overview_xsection.png")

uc_overview = quickplot(obs, ihacres_node.outflow, climate, "Split V NmKGE", true)

plot(bl_overview, uc_overview, layout=(2, 1), size=(850, 700))
savefig("overview_qq.png")

# Streamfall.temporal_cross_section(climate_data.Date, h_data, n_data; period=monthday)
# savefig("temporal_xsection_monthday_ME.png")

# Streamfall.temporal_cross_section(climate_data.Date, h_data, n_data; period=month, func=Streamfall.mKGE)
# savefig("test_mKGE.png")

# Displaying results and saving figure
# plot(h_data,
#      legend=:bottomleft,
#      title="Calibrated IHACRES\n(RMSE: $(rmse); NSE: $(nse))",
#      label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

# display(plot!(n_data, label="IHACRES"))

# savefig("calibrated_example.png")
