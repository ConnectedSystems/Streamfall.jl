using DataFrames, CSV
using Statistics
using BlackBoxOptim

using Infiltrator
using ModelParameters


include("./network_creation.jl")

climate_data = DataFrame!(CSV.File("../test/data/climate/climate_historic.csv", 
                          comment="#",
                          dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File("../test/data/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File("../test/data/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

# Subset to same range
last_date = hist_dam_levels.Date[end]
climate_data = climate_data[climate_data.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[hist_dam_releases.Date .<= last_date, :]

climate = Climate(climate_data, "_rain", "_evap")


function obj_func(params)
    global climate
    global hist_dam_levels
    global hist_dam_releases
    global mg
    global g
    global v_id

    node = get_prop(mg, v_id, :node)
    node = deepcopy(node)
    update_params!(node, params)  # node = ModelParameters.update(node, params)
    set_prop!(mg, v_id, :node, node)

    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        run_node!(mg, g, v_id, climate, ts; water_order=hist_dam_releases)
    end

    # Calculate score (NSE)
    hist_levels = hist_dam_levels[:, "Dam Level [mAHD]"]
    score = -(1.0 - sum((node.level .- hist_levels).^2) / sum((hist_levels .- mean(hist_levels)).^2))

    # reset to clear stored values
    # reset!(get_prop(mg, v_id, :node))

    return score
end


function calibrate(g, mg, v_id)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(g, mg, ins)
        end
    end

    @infiltrate

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, mg, g, v_id, discharge, levels)

    target_node = Model(get_prop(mg, v_id, :node))
    score = bboptimize(opt_func; SearchRange=collect(target_node.bounds))

    return score
end


match = collect(filter_vertices(mg, :name, "406000"))
v_id = match[1]
# target_node = Model(get_prop(mg, v_id, :node))

@info calibrate(g, mg, v_id)
