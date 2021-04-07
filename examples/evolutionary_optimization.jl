using Distributed
using Evolutionary

using Infiltrator

# addprocs(2, exeflags="--project=../")

@everywhere begin
    # Ensure dependent data and packages are available
    using DataFrames, CSV
    using Statistics
    using Evolutionary

    using ModelParameters
    using LightGraphs, MetaGraphs
    using YAML
    using Streamfall


    network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
    g, mg = create_network("Example Network", network)
    inlets, outlets = find_inlets_and_outlets(g)

    # @info "Network has the following inlets and outlets:" inlets outlets

    climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                              comment="#",
                              dateformat="YYYY-mm-dd"))

    hist_dam_levels = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
    hist_dam_releases = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

    inlet_levels = DataFrame!(CSV.File("../test/data/campaspe/gauges/406219_edited.csv", dateformat="YYYY-mm-dd"))
    inlet_flows = DataFrame!(CSV.File("../test/data/campaspe/gauges/406219_outflow_edited.csv", dateformat="YYYY-mm-dd"))

    # Subset to same range
    first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1], inlet_flows.Date[1])
    last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end], inlet_flows.Date[end])

    # @info "Date ranges:" first_date last_date

    climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
    hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]
    hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]
    # inlet_levels = inlet_levels[first_date .<= inlet_levels.Date .<= last_date, :]
    inlet_flows = inlet_flows[first_date .<= inlet_flows.Date .<= last_date, :]

    hist_data = Dict(
        "406000" => hist_dam_levels[:, "Dam Level [mAHD]"],
        "406219" => inlet_flows[:, "406219_outflow_[ML]"]
    )

    climate = Climate(climate_data, "_rain", "_evap")

    function obj_func(params, climate, mg, g, v_id, obs_data)

        node = get_prop(mg, v_id, :node)
        update_params!(node, params...)
    
        timesteps = sim_length(climate)
        for ts in (1:timesteps)
            run_node!(mg, g, v_id, climate, ts; water_order=hist_dam_releases)
        end

        if node.node_id == "406000"
            node_data = node.level[10:end]
        else
            node_data = node.outflow[10:end]
        end

        h_data = obs_data[10:end]
    
        # Calculate score (NSE)
        NSE = 1 - sum((h_data .- node_data).^2) / sum((h_data .- mean(h_data)).^2)
    
        # Normalized NSE so that score ranges from 0 to 1. NNSE of 0.5 is equivalent to NSE = 0.
        NNSE = 1 / (2 - NSE)
    
        # Swap signs as we want to minimize
        score = -NNSE
    
        # reset to clear stored values
        reset!(node)
    
        # Borg method expects tuple to be returned
        # return (score, )
        return score
    end
end


function calibrate(g, mg, v_id, climate, calib_data)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(g, mg, ins, climate, calib_data)
        end
    end

    node = get_prop(mg, v_id, :node)
    node_id = node.node_id

    obs_data = calib_data[node_id]

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, mg, g, v_id, obs_data)

    # Get all parameters
    x0, param_bounds = param_info(node)
    lower, upper = collect(zip(param_bounds...))
    lower, upper = collect(lower), collect(upper)

    cnst = BoxConstraints(lower, upper)
    res = Evolutionary.optimize(opt_func, cnst, x0, 
                                CMAES() # DE()
    )

    bs = Evolutionary.minimizer(res)
    @info "Calibrated $(v_id) ($(node_id)), with score: $(Evolutionary.minimum(res))"
    @info "Best Params:" bs
    # @info converged(res)

    # Update node with calibrated parameters
    update_params!(get_prop(mg, v_id, :node), bs...)

    return res
end


match = collect(filter_vertices(mg, :name, "406000"))
v_id = match[1]

@info "Starting calibration..."
res = calibrate(g, mg, v_id, climate, hist_data)

@info "Score:" Evolutionary.minimum(res)
@info "Best Params:" Evolutionary.minimizer(res)

node = get_prop(mg, v_id, :node)
@info node


using Plots

timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(mg, g, 1, climate, ts)  # ; water_order=hist_dam_releases
end

match = collect(filter_vertices(mg, :name, "406219"))
inlet_id = match[1]

node = get_prop(mg, inlet_id, :node)
h_data = inlet_flows[:, "406219_outflow_[ML]"]
n_data = node.outflow

plot(n_data)
plot!(h_data)

NSE = 1.0 - sum((h_data .- n_data).^2) / sum((h_data .- mean(h_data)).^2)
NNSE = 1.0 / (2.0 - NSE)

@info "NNSE:" NNSE


# timesteps = sim_length(climate)
# for ts in (1:timesteps)
#     run_node!(mg, g, 2, climate, ts; water_order=hist_dam_releases)
# end

# node = get_prop(mg, 2, :node)
# plot(node.level)
# plot!(hist_dam_levels[:, "Dam Level [mAHD]"])
