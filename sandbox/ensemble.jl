using Streamfall
using Statistics, CSV, DataFrames
using Distributed, BlackBoxOptim
using Base.Iterators
using Plots
import Streamfall.Analysis: TemporalCrossSection


# Spin up workers if needed
if nprocs() == 1
    addprocs(4, exeflags="--project=..")
    @everywhere begin
        # Activate environment and precompile
        # using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
        # Pkg.instantiate(); Pkg.precompile()

        using Streamfall
        using Statistics, CSV, DataFrames
    end
end


@everywhere begin

function v_metric(obs, sim; sp=365, metric=Streamfall.NmKGE)
    res = Streamfall.naive_split_metric(obs, sim, sp, metric)
    m = mean(res)
    cv = std(res) / m

    V = (m*0.9) + ((ℯ^-cv)*0.1)  #  (ℯ .- 1.5)

    return V
end


function v_mod_metric(obs, sim; sp=365, metric=Streamfall.NmKGE)
    res = Streamfall.naive_split_metric(obs, sim, sp, metric)
    m = mean(res)
    sd = std(res)

    V = (m*0.9) + ((ℯ^-sd)*0.1)

    return V
end


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

# Create nodes/ensemble
ihacres_node = create_node(BilinearNode, "410730", 129.2)
gr4j_node = create_node(GR4JNode, "410730", 129.2)
# ensemble_node = BaseEnsemble([ihacres_node, gr4j_node], weights=[0.5, 0.5])


burn_in = 365

split_v = (obs, sim) -> 1.0 - v_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)
split_mod_v = (obs, sim) -> 1.0 - v_mod_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)

end  # everywhere block


name_nodes = zip(["IHACRES", "GR4J"], [ihacres_node, gr4j_node])
name_metrics = zip(["split_v", "split_mod_v", "NnpKGE", "NmKGE"], [split_v, split_mod_v, Streamfall.NnpKGE, Streamfall.NmKGE])
node_metric = [NamedTuple{(:name, :node, :met_name, :metric)}([name, node, met_name, metric]) for ((name, node), (met_name, metric)) in product(name_nodes, name_metrics)]


# res, opt = calibrate!(ihacres_node, climate, hist_streamflow; metric=split_v, MaxTime=1200.0, PopulationSize=500)
# res, opt = calibrate!(gr4j_node, climate, hist_streamflow; metric=split_v, MaxTime=1200.0, PopulationSize=500)


# p_res = pmap((n, m) -> calibrate!(n, climate, hist_streamflow; metric=m, MaxTime=1200.0))

p_res = pmap((x) -> calibrate!(x.node, climate, hist_streamflow; metric=x.metric, MaxTime=600.0, TraceMode=:silent), node_metric)


function postprocess(x, ω)
    x = copy(x)
    idx = (x .+ ω) .>= 0.0
    x[idx] = x[idx] .+ ω[idx]

    # idx2 = (x .+ ω) .< 0.0
    # x[idx2] .= 0.0
    return x
end


c_date = climate_data.Date
obs = hist_streamflow[:, "410730"]
figs = []
outflows = []
display_timeframe = 13330:(13330+(365*1))
for (x, (r, o)) in zip(node_metric, p_res)
    update_params!(x.node, best_candidate(r)...)
    reset!(x.node)
    run_node!(x.node, climate)

    t_name = x.name
    m_name = x.met_name
    push!(outflows, x.node.outflow)

    mae = round(Streamfall.MAE(obs, x.node.outflow), digits=4)
    push!(figs, Streamfall.temporal_cross_section(c_date, obs, x.node.outflow; yscale=:log10, title="$t_name ($m_name | $mae)"))

    mod_xsect_res = TemporalCrossSection(c_date, obs, x.node.outflow)
    offset_vals = mod_xsect_res.offsets

    pp_outflow = postprocess(x.node.outflow, offset_vals)
    mae = round(Streamfall.MAE(obs, pp_outflow), digits=4)
    push!(figs, Streamfall.temporal_cross_section(c_date, obs, pp_outflow; yscale=:log10, title="Post-Processed | $mae"))

    t_plot = plot(c_date[display_timeframe], obs[display_timeframe], label="Observed")
    plot!(c_date[display_timeframe], x.node.outflow[display_timeframe], label="Modeled")

    push!(figs, t_plot)

    pp_t_plot = plot(c_date[display_timeframe], obs[display_timeframe], label="Observed")
    plot!(c_date[display_timeframe], pp_outflow[display_timeframe], label="Post-Processed")

    push!(figs, pp_t_plot)
end

# ensemble_node = BaseEnsemble([ihacres_node, gr4j_node], weights=[0.5, 0.5])
# calibrate!(ensemble_node, climate, hist_streamflow; metric=metric, MaxTime=120.0, TraceMode=:silent)


plot(figs..., size=((300*4),(250*8)), layout=(8,4), dpi=300)  # , link=:y
savefig("metric_calibration_comparison_2022-04-31__90_10_weight_1_year.png")

# TODO: Try heatmap of mean error over day of year

# TODO: Reincorporate ensemble

# TODO: Try out state-based implementation

# quickplot(obs, outflows[2], climate, "", false; burn_in=13330, limit=13330+3650)

# # p_res will have pairs of res and opt, where res are the results and opt are the options used
# # bs = best_candidate(res)  # extract best parameters
# # update_params!(node, bs...)  # update node with parameters

# # calibrate!(ensemble_node, climate, hist_streamflow; metric=split_v, MaxTime=1500.0, PopulationSize=500)

# run_node!(ensemble_node, climate)

# obs = hist_streamflow[:, "410730"]

# quickplot(obs, ensemble_node, climate; burn_in=365, metric=split_v)
# quickplot(obs, ihacres_node, climate; burn_in=365, metric=split_v)

# c_date = climate_data.Date

# #########
# mod_xsect_res = TemporalCrossSection(c_date, obs, ensemble_node.outflow)
# offset_vals = mod_xsect_res.offsets

# # fig2 = Streamfall.temporal_cross_section(c_date, obs, ensemble_node.outflow; yscale=:log10, title="NnpKGE")

# function postprocess(x, ω)
#     x = copy(x)
#     idx = (x .+ ω) .>= 0.0
#     x[idx] = x[idx] .+ ω[idx]

#     # idx2 = (x .+ ω) .< 0.0
#     # x[idx2] .= 0.0
#     return x
# end

# #########


# fig1 = Streamfall.temporal_cross_section(c_date, obs, ihacres_node.outflow; yscale=:log10, title="IHACRES")
# fig2 = Streamfall.temporal_cross_section(c_date, obs, gr4j_node.outflow; yscale=:log10, title="GR4J")
# fig3 = Streamfall.temporal_cross_section(c_date, obs, ensemble_node.outflow; yscale=:log10, title="Ensemble (IHACRES-GR4J)")
# fig4 = Streamfall.temporal_cross_section(c_date, obs, postprocess(ensemble_node.outflow, offset_vals); yscale=:log10, title="Post-Processed (IHACRES-GR4J)")

# plot(fig1, fig2, fig3, fig4, size=(1000,800), layout=(2,2), link=:y, dpi=300)

# savefig("ensemble_pp_popsize_500_maxtime_900_modded_v_metric.png")

# quickplot(obs, postprocess(ensemble_node.outflow, offset_vals), c_date)
# savefig("ensemble_qp_pp_popsize_500_modded_v_metric.png")

# quickplot(obs, ensemble_node.outflow, c_date)
# savefig("ensemble_qp_ensemble_popsize_500_modded_v_metric.png")

# quickplot(obs, ihacres_node.outflow, c_date)
# savefig("ensemble_qp_ihacres_popsize_500_modded_v_metric.png")

# quickplot(obs, gr4j_node.outflow, c_date)
# savefig("ensemble_qp_gr4j_popsize_500_modded_v_metric.png")