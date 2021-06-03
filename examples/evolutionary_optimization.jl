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
    sn = create_network("Example Network", network)
    inlets, outlets = find_inlets_and_outlets(sn)

    # @info "Network has the following inlets and outlets:" inlets outlets

    climate_data = DataFrame!(CSV.File("../test/data/campaspe/climate/climate_historic.csv", 
                              comment="#",
                              dateformat="YYYY-mm-dd"))

    hist_dam_levels = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_levels_for_fit.csv", dateformat="YYYY-mm-dd"))
    hist_dam_releases = DataFrame!(CSV.File("../test/data/campaspe/dam/historic_releases.csv", dateformat="YYYY-mm-dd"))

    inlet_levels = DataFrame!(CSV.File("../test/data/campaspe/gauges/406219_edited.csv", dateformat="YYYY-mm-dd"))

    # Subset to same range
    first_date = max(hist_dam_levels.Date[1], hist_dam_releases.Date[1])
    last_date = min(hist_dam_levels.Date[end], hist_dam_releases.Date[end])

    climate_data = climate_data[first_date .<= climate_data.Date .<= last_date, :]
    hist_dam_releases = hist_dam_releases[first_date .<= hist_dam_releases.Date .<= last_date, :]
    hist_dam_levels = hist_dam_levels[first_date .<= hist_dam_levels.Date .<= last_date, :]

    hist_data = Dict(
        "406000" => hist_dam_levels[:, "Dam Level [mAHD]"]
    )

    climate = Climate(climate_data, "_rain", "_evap")

    function obj_func(params, climate, sn, v_id, next_vid, calib_data)

        this_node = get_prop(sn, v_id, :node)
        update_params!(this_node, params...)

        next_node = get_prop(sn, next_vid, :node)
    
        timesteps = sim_length(climate)
        for ts in (1:timesteps)
            run_node!(sn, next_vid, climate, ts; water_order=hist_dam_releases)
        end

        if next_node.node_id == "406000"
            node_data = next_node.level
            h_data = calib_data[next_node.node_id]
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
        reset!(sn)
    
        # Borg method expects tuple to be returned
        # return (score, )
        return score
    end
end


function calibrate(sn, v_id, climate, calib_data)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(sn, ins, climate, calib_data)
        end
    end

    outs = outneighbors(sn, v_id)
    @assert length(outs) == 1 || throw("Streamfall currently only supports a single outlet.")
    outs = outs[1]

    this_node = get_prop(sn, v_id, :node)
    next_node = get_prop(sn, outs, :node)
    next_id = next_node.node_id

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, sn, v_id, outs, calib_data)

    # Get node parameters
    x0, param_bounds = param_info(this_node; with_level=false)
    lower, upper = collect(zip(param_bounds...))
    lower, upper = collect(lower), collect(upper)

    cnst = BoxConstraints(lower, upper)
    res = Evolutionary.optimize(opt_func, cnst, x0, 
                                DE(), # CMAES()
                                Evolutionary.Options(iterations=Int(1e5), abstol=1e-32)
    )

    bs = Evolutionary.minimizer(res)
    @info "Calibrated $(v_id) ($(this_node.node_id)), with score: $(Evolutionary.minimum(res))"
    @info "Best Params:" bs

    # Update node with calibrated parameters
    update_params!(get_prop(sn, v_id, :node), bs...)

    return res
end


match = collect(filter_vertices(sn, :name, "406219"))
v_id = match[1]

@info "Starting calibration..."
res = calibrate(sn, v_id, climate, hist_data)

@info best_fitness(res)
@info best_candidate(res)

node = get_prop(sn, v_id, :node)
@info node


using Plots

match = collect(filter_vertices(sn, :name, "406000"))
dam_id = match[1]

timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(sn, dam_id, climate, ts; water_order=hist_dam_releases)
end

dam_node = get_prop(sn, dam_id, :node)
h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

@info "NNSE:" Streamfall.NNSE(h_data, n_data)
@info "RMSE:" Streamfall.RMSE(h_data, n_data)

plot(h_data)
plot!(n_data)
