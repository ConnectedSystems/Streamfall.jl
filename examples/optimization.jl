using Distributed
using BlackBoxOptim

# addprocs(2, exeflags="--project=../")

@everywhere begin
    # Ensure dependent data and packages are available
    using BlackBoxOptim
    using DataFrames, CSV
    using Statistics

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


    function obj_func(params, climate, mg, g, v_id, next_vid, calib_data)

        this_node = get_prop(mg, v_id, :node)
        update_params!(this_node, params...)

        next_node = get_prop(mg, next_vid, :node)
    
        timesteps = sim_length(climate)
        for ts in (1:timesteps)
            run_node!(mg, g, next_vid, climate, ts; water_order=hist_dam_releases)
        end

        if next_node.node_id == "406000"
            node_data = next_node.level
            h_data = calib_data[next_node.node_id]
        else
            node_data = next_node.outflow
            h_data = calib_data[this_node.node_id]
        end
    
        # Calculate score (NSE)
        NNSE = Streamfall.NNSE(h_data, node_data)
    
        # Swap signs as we want to minimize
        score = -NNSE
    
        # RMSE = (sum((node_levels .- h_levels).^2)/length(node_levels))^0.5
        # score = RMSE
    
        # reset to clear stored values
        reset!(this_node)
    
        # Borg method expects tuple to be returned
        # return (score, )
        return score
    end
end


function calibrate(mg, g, v_id, climate, calib_data)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(mg, g, ins, climate, calib_data)
        end
    end

    outs = outneighbors(g, v_id)
    @assert length(outs) == 1 || throw("Streamfall currently only supports a single outlet. ($(length(outs)))")
    outs = outs[1]

    this_node = get_prop(mg, v_id, :node)
    next_node = get_prop(mg, outs, :node)
    next_id = next_node.node_id

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, mg, g, v_id, outs, calib_data)

    # Get node parameters
    x0, param_bounds = param_info(this_node; with_level=false)

    # opt = bbsetup(opt_func; SearchRange=param_bounds,
    #               Method=:borg_moea,
    #               ϵ=0.05,
    #               FitnessScheme=ParetoFitnessScheme{1}(is_minimizing=true),
    #               MaxTime=600.0,  #spend 10 minutes calibrating each node
    #               TraceMode=:silent,
    #               Workers=workers())
    # NThreads=Threads.nthreads()-1
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  Method=:adaptive_de_rand_1_bin_radiuslimited,
                  # MaxTime=600.0,  #spend 10 minutes calibrating each node
                  TraceMode=:silent,
                  PopulationSize=100,
                  ftol=-1e-10,
                  Workers=workers())
    
    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(this_node.node_id)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(get_prop(mg, v_id, :node), bs...)

    return res
end


match = collect(filter_vertices(mg, :name, "406219"))
v_id = match[1]

@info "Starting calibration..."
res = calibrate(mg, g, v_id, climate, hist_data)

@info best_fitness(res)
@info best_candidate(res)

node = get_prop(mg, v_id, :node)
@info node


using Plots


match = collect(filter_vertices(mg, :name, "406000"))
dam_id = match[1]

timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(mg, g, dam_id, climate, ts; water_order=hist_dam_releases)
end

dam_node = get_prop(mg, dam_id, :node)
h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

@info "NNSE:" Streamfall.NNSE(h_data, n_data)
@info "RMSE:" Streamfall.RMSE(h_data, n_data)


plot(h_data)
plot!(n_data)



# timesteps = sim_length(climate)
# for ts in (1:timesteps)
#     run_node!(mg, g, 2, climate, ts; water_order=hist_dam_releases)
# end

# node = get_prop(mg, 2, :node)
# plot(node.level)
# plot!(hist_dam_levels[:, "Dam Level [mAHD]"])
