"""Another calibration example, this time using split-mKGE.
"""

# Import common packages and functions
include("_obj_func_definition.jl")


"""Another Example calibration function.

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

    # Need to get outlet if performance of current node is dependent on next node.
    # e.g., if you want to calibrate against downstream dam levels.
    # There should only be one outlet
    out_node_id = outlets(sn, v_id)[1]
    this_node = sn[v_id]

    # Create new optimization function (see definition inside `_obj_func_definition.jl`)
    # opt_func = x -> obj_func(x, climate, sn, v_id, outs, calib_data)
    opt_func = x -> alt_obj_func(x, climate, sn, v_id, out_node_id, calib_data)

    # Get node parameters (default values and bounds)
    p_names, x0, param_bounds = param_info(this_node; with_level=false)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  Method=:adaptive_de_rand_1_bin_radiuslimited,
                  MaxTime=2400.0,  # time in seconds to spend
                  TraceInterval=30.0,
                  PopulationSize=75,
                  )

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(this_node.node_id)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(this_node, bs...)

    return res, opt
end


v_id, node = sn["406219"]
@info "Starting calibration..."
res, opt = calibrate(sn, v_id, climate, hist_data)

best_params = best_candidate(res)

@info best_fitness(res)
@info best_params

reset!(node)

using Plots

update_params!(node, best_params...)

dam_id, dam_node = sn["406000"]
run_node!(sn, dam_id, climate; extraction=hist_dam_releases)

h_data = hist_data["406000"]
n_data = dam_node.level

split_mKGE = Streamfall.naive_split_metric(h_data, n_data; n_members=365, metric=Streamfall.NmKGE, comb_method=mean)
rmse = Streamfall.RMSE(h_data, n_data)

split_mKGE = round(split_mKGE, digits=4)
rmse = round(rmse, digits=4)

plot(h_data,
     legend=:bottomleft,
     title="Calibrated IHACRES\n(Split NmKGE: $(split_mKGE); RMSE: $(rmse))",
     label="Historic", xlabel="Day", ylabel="Dam Level [mAHD]")

plot!(n_data, label="IHACRES")

# RMSE calibrated best candidate: [124.069, 1.97, 0.823727, 1.95943, 8.1387, 0.0988389, 9.48091, 0.696038]
# Split RMSE calibrated best candidate
savefig("calibrated_Split_NmKGE.png")