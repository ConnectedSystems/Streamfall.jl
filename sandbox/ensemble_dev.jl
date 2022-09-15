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

        include("IHACRESCMDStateNode.jl")
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

state_node = create_node(CMDStateNode, "410730", 129.2)


burn_in = 1826

split_v = (obs, sim) -> 1.0 - v_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)
split_mod_v = (obs, sim) -> 1.0 - v_mod_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)

end  # everywhere block


name_nodes = zip(["IHACRES", "GR4J", "cmd_state"], [ihacres_node, gr4j_node, state_node])
name_metrics = zip(["split_v", "split_mod_v", "NnpKGE", "NmKGE"], [split_v, split_mod_v, Streamfall.NnpKGE, Streamfall.NmKGE])
node_metric = [NamedTuple{(:name, :node, :met_name, :metric)}([name, node, met_name, metric]) for ((name, node), (met_name, metric)) in product(name_nodes, name_metrics)]

p_res = pmap((x) -> calibrate!(x.node, climate, hist_streamflow; metric=x.metric, MaxTime=1800.0, TraceMode=:silent), node_metric)


function postprocess(x, ω)
    x = copy(x)
    idx = (x .+ ω) .>= 0.0
    x[idx] = x[idx] .+ ω[idx]

    # x .= max.(x .+ ω, 0.0)
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


# split_mod_v2 = (obs, sim) -> 1.0 - v_mod_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)


## Ensemble modelling section
# ensemble_node = BaseEnsemble([node_metric[3].node, node_metric[4].node], weights=[0.5, 0.5])
# reset!(ensemble_node)
# e_res, e_opt = calibrate!(ensemble_node, climate, hist_streamflow; metric=split_mod_v, MaxTime=1200.0, TraceMode=:silent)

# update_params!(ensemble_node, best_candidate(e_res)...)
# reset!(ensemble_node)
# run_node!(ensemble_node, climate)

# push!(outflows, ensemble_node.outflow)

# mae = round(Streamfall.MAE(obs, ensemble_node.outflow), digits=4)
# push!(figs, Streamfall.temporal_cross_section(c_date, obs, ensemble_node.outflow; yscale=:log10, title="Ensemble (split_mod_v | $mae)"))

# mod_xsect_res = TemporalCrossSection(c_date, obs, ensemble_node.outflow)
# offset_vals = mod_xsect_res.offsets

# pp_outflow = postprocess(ensemble_node.outflow, offset_vals)
# mae = round(Streamfall.MAE(obs, pp_outflow), digits=4)
# push!(figs, Streamfall.temporal_cross_section(c_date, obs, pp_outflow; yscale=:log10, title="Post-Processed | $mae"))

# t_plot = plot(c_date[display_timeframe], obs[display_timeframe], label="Observed")
# plot!(c_date[display_timeframe], ensemble_node.outflow[display_timeframe], label="Modeled")
# push!(figs, t_plot)

# pp_t_plot = plot(c_date[display_timeframe], obs[display_timeframe], label="Observed")
# plot!(c_date[display_timeframe], pp_outflow[display_timeframe], label="Post-Processed")
# push!(figs, pp_t_plot)


plot(figs..., size=((350*4),(300*12)), layout=(12,4), dpi=300)  # , link=:y
savefig("figs/metric_calibration_comparison_2022-05-01b__90_10_weight_1_year.png")




# TODO: Try heatmap of mean error over day of year

# TODO: Reincorporate ensemble [done]

# TODO: Try out state-based implementation [done]

# TODO: Try out state-based implementation with borg moea, targeting each catchment state
