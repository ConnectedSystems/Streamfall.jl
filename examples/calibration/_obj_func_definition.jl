DATA_PATH = "../../test/data/campaspe/"

# Ensure dependent data and packages are available
using DataFrames, CSV
using Statistics
using Evolutionary

using ModelParameters
using LightGraphs, MetaGraphs
using YAML
using Streamfall


network = YAML.load_file(joinpath(DATA_PATH, "campaspe_network.yml"))
mg, g = create_network("Example Network", network)
inlets, outlets = find_inlets_and_outlets(g)

# @info "Network has the following inlets and outlets:" inlets outlets

climate_data = DataFrame!(CSV.File(joinpath(DATA_PATH, "climate/climate_historic.csv"),
                            comment="#",
                            dateformat="YYYY-mm-dd"))

hist_dam_levels = DataFrame!(CSV.File(joinpath(DATA_PATH, "dam/historic_levels_for_fit.csv"), dateformat="YYYY-mm-dd"))
hist_dam_releases = DataFrame!(CSV.File(joinpath(DATA_PATH, "dam/historic_releases.csv"), dateformat="YYYY-mm-dd"))

# Subset to same range
first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1], climate_data.Date[1])
last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end], climate_data.Date[end])

climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]
hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]

hist_data = Dict(
    "406000" => hist_dam_levels[:, "Dam Level [mAHD]"]
)

climate = Climate(climate_data, "_rain", "_evap")


function obj_func(params, climate, mg, g, v_id, next_vid, calib_data)

    this_node = get_node(mg, v_id)
    update_params!(this_node, params...)

    next_node = get_node(mg, next_vid)

    timesteps = sim_length(climate)
    for ts in (1:timesteps)
        Streamfall.run_node!(mg, g, next_vid, climate, ts; water_order=hist_dam_releases)
    end

    if next_node.node_id == "406000"
        node_data = next_node.level
        h_data = calib_data[next_node.node_id]
    elseif this_node.node_id == "406000"
        node_data = this_node.level
        h_data = calib_data[this_node.node_id]
    else
        node_data = this_node.outflow
        h_data = calib_data[this_node.node_id]
    end

    # Calculate score (NNSE; 0 to 1)
    NNSE = Streamfall.NNSE(h_data, node_data)

    # Switch fitness direction as we want to minimize
    score = 1.0 - NNSE

    # RMSE = Streamfall.RMSE(h_data, node_data)
    # score = RMSE

    # reset to clear stored values
    reset!(mg, g)

    # Borg method expects tuple to be returned
    # return (score, )
    return score
end