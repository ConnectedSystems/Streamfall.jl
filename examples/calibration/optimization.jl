using Distributed
using BlackBoxOptim

# if nworkers() < 3
#     addprocs(2, exeflags="--project=../")
# end

include("_obj_func_definition.jl")


function calibrate(mg, g, v_id, climate, calib_data)

    inlets = inneighbors(g, v_id)
    if !isempty(inlets)
        for ins in inlets
            calibrate(mg, g, ins, climate, calib_data)
        end
    end

    outs = outneighbors(g, v_id)
    @assert length(outs) <= 1 || throw("Streamfall currently only supports a single outlet. ($(length(outs)))")
    outs = outs[1]

    this_node = get_node(mg, v_id)

    # Create new optimization function
    opt_func = x -> obj_func(x, climate, mg, g, v_id, outs, calib_data)

    # Get node parameters
    x0, param_bounds = param_info(this_node; with_level=false)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  Method=:adaptive_de_rand_1_bin_radiuslimited,
                  MaxTime=3600.0,  # time in seconds to spend
                  TraceInterval=10.0,
                  PopulationSize=100,
                  # Workers=workers()
                  )
    
    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(this_node.node_id)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(this_node, bs...)

    return res, opt
end


v_id, node = get_gauge(mg, "406219")
@info "Starting calibration..."
res, opt = calibrate(mg, g, v_id, climate, hist_data)

@info best_fitness(res)
@info best_candidate(res)

@info node


using Plots

dam_id, dam_node = get_gauge(mg, "406000")
timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(mg, g, dam_id, climate, ts; water_order=hist_dam_releases)
end

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
