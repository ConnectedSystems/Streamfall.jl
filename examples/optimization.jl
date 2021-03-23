using DataFrames, CSV
using Statistics
using BlackBoxOptim

using Infiltrator
using ModelParameters
using LightGraphs, MetaGraphs


include("./network_creation.jl")

climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

inlet_levels = DataFrame!(CSV.File("../test/data/campaspe/gauges/406219_edited.csv", dateformat="YYYY-mm-dd"))


# Subset to same range
first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1], inlet_levels.Date[1])
last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end], inlet_levels.Date[end])

@info "Date ranges:" first_date last_date

climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]
hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]
inlet_levels = inlet_levels[first_date .<= inlet_levels.Date .<= last_date, :]

hist_levels = Dict(
    "406000" => hist_dam_levels,
    "406219" => inlet_levels
)

climate = Climate(climate_data, "_rain", "_evap")


function obj_func(params, climate, mg, g, v_id, hist_levels)
    global hist_dam_releases

    orig_node = get_prop(mg, v_id, :node)
    node = deepcopy(orig_node)
    update_params!(node, params...)
    set_prop!(mg, v_id, :node, node)

    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(mg, g, v_id, climate, ts; water_order=hist_dam_releases)
    end

    node_id = node.node_id
    obs_levels = hist_levels[node_id]
    if node_id == "406000"
        h_levels = obs_levels[:, "Dam Level [mAHD]"]
    else
        h_levels = obs_levels[:, "$(node_id)_level"]
    end

    node_levels = node.level[1000:end]
    h_levels = h_levels[1000:end]

    # Calculate score (NSE)
    NSE = 1 - sum((h_levels .- node_levels).^2) / sum((h_levels .- mean(h_levels)).^2)

    # Normalized NSE so that score ranges from 0 to 1. NNSE of 0.5 is equivalent to NSE = 0.
    NNSE = 1 / (2 - NSE)

    # Swap signs as we want to minimize
    score = -NNSE

    # RMSE = (sum((node_levels .- h_levels).^2)/length(node_levels))^0.5
    # score = RMSE

    # reset to clear stored values
    # set_prop!(mg, v_id, :node, orig_node)
    reset!(node)

    # Borg method expects tuple to be returned
    return (score, )
end


function calibrate(g, mg, v_id)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(g, mg, ins)
        end
    end

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, mg, g, v_id, hist_levels)

    node = get_prop(mg, v_id, :node)
    target_node = Model(node)

    # Add area as static parameter (assume known)
    # param_bounds = nothing
    # try
    #     param_bounds = vcat((node.area, node.area), target_node.bounds..., [x.bounds for x in node.level_params])
    # catch
    #     param_bounds = collect(target_node.bounds)
    # end
    param_bounds = collect(target_node.bounds)

    res = bboptimize(opt_func; SearchRange=param_bounds,
                       Method=:borg_moea,
                       FitnessScheme=ParetoFitnessScheme{1}(is_minimizing=true),
                       MaxTime=180.0,
                       TraceMode = :silent)

    @info "Calibrated $(v_id), with score: $(best_fitness(res))"

    bs = best_candidate(res)
    update_params!(get_prop(mg, v_id, :node), bs...)

    return res
end


match = collect(filter_vertices(mg, :name, "406000"))
v_id = match[1]
# target_node = Model(get_prop(mg, v_id, :node))

res = calibrate(g, mg, v_id)

@info best_fitness(res)
@info best_candidate(res)

node = get_prop(mg, v_id, :node)
@info node


using Plots

timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(mg, g, 1, climate, ts; water_order=hist_dam_releases)
end

node = get_prop(mg, 1, :node)
h_level = inlet_levels[:, "406219_level"]
plot(node.level)
plot!(h_level)

NSE = 1 - sum((h_level .- node.level).^2) / sum((h_level .- mean(h_level)).^2)

@info NSE


timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(mg, g, 2, climate, ts; water_order=hist_dam_releases)
end

node = get_prop(mg, 2, :node)
plot(node.level)
plot!(hist_dam_levels[:, "Dam Level [mAHD]"])
