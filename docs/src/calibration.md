# Example calibration


```julia
# Import common packages and functions
# This is shown under Calibration Setup
include("_obj_func_definition.jl")


"""Example calibration function.

Illustrate model calibration using the BlackBoxOptim package.
"""
function calibrate(sn, v_id, climate, calib_data)

    ins = inlets(sn, v_id)

    # Recurse through and calibrate all nodes upstream
    if !isempty(ins)
        for nid in ins
            calibrate(sn, nid, climate, calib_data)
        end
    end

    this_node = get_node(sn, v_id)

    # Create new optimization function (see definition inside Calibration Setup)
    opt_func = x -> obj_func(x, climate, sn, v_id, calib_data)

    # Get node parameters (default values and bounds)
    x0, param_bounds = param_info(this_node; with_level=false)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  Method=:adaptive_de_rand_1_bin_radiuslimited,
                  MaxTime=300.0,  # time in seconds to spend
                  TraceInterval=30.0,
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


v_id, node = get_gauge(sn, "406219")
@info "Starting calibration..."
res, opt = calibrate(sn, v_id, climate, hist_data)

best_params = best_candidate(res)

@info best_fitness(res)
@info best_params


update_params!(node, best_params...)
Streamfall.run_node!(sn, v_id, climate)

h_data = hist_data["406219"]
n_data = node.outflow

@info "Outflow NNSE:" Streamfall.NNSE(h_data, n_data)
@info "Outflow RMSE:" Streamfall.RMSE(h_data, n_data)
reset!(node)

dam_id, dam_node = get_gauge(sn, "406000")
timesteps = sim_length(climate)
for ts in (1:timesteps)
    run_node!(sn, dam_id, climate, ts; water_order=hist_dam_releases)
end

h_data = hist_dam_levels[:, "Dam Level [mAHD]"]
n_data = dam_node.level

@info "Downstream Dam Level NNSE:" Streamfall.NNSE(h_data, n_data)
@info "Downstream Dam Level RMSE:" Streamfall.RMSE(h_data, n_data)


# Best candidate found: [54.9098, 0.135862, 1.22086, 2.99995, 0.309896, 0.0861618, 0.977643, 0.869782]
```