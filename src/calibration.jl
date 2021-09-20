# Ensure dependent data and packages are available
using Distributed, BlackBoxOptim, Serialization


"""Calibrate specified node in network."""
function obj_func(params, climate::Climate, sn::StreamfallNetwork, v_id::Int, calib_data::Array; 
                  metric::Function, inflow=nothing, extraction=nothing, exchange=nothing)
    return obj_func(params, climate, sn[v_id], calib_data; metric=metric, inflow=inflow, extraction=extraction, exchange=exchange)
end


"""Calibrate current node."""
function obj_func(params, climate::Climate, node::NetworkNode, calib_data::Array; 
                  metric::Function, inflow=nothing, extraction=nothing, exchange=nothing)
    update_params!(node, params...)

    run_node!(node, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    n_data = node.outflow
    score = metric(calib_data, n_data)

    # Reset to clear stored values
    reset!(node)

    return score
end


"""
    calibrate!(sn, v_id, climate, calib_data,
               metric::Function=Streamfall.RMSE;
               kwargs...)

Calibrate a given node, recursing upstream, using the BlackBoxOptim package.

# Arguments
- `sn::StreamfallNetwork` : Network
- `v_id::Int` : node identifier
- `climate::Climate` : Climate data
- `calib_data::Array` : calibration data for target node by its id
- `metric::Function` : Optimization function to use. Defaults to RMSE.
"""
function calibrate!(sn::StreamfallNetwork, v_id::Int64, climate::Climate, calib_data::DataFrame;
                    metric::Function=Streamfall.RMSE,
                    kwargs...)
    calibrate!(sn, v_id, climate, calib_data[:, sn[v_id].name]; metric=metric, kwargs...)
end


"""
    calibrate!(sn, v_id, climate, calib_data,
               metric::Function=Streamfall.RMSE;
               kwargs...)

Calibrate a given node, recursing upstream, using the BlackBoxOptim package.

# Arguments
- `sn::StreamfallNetwork` : Network
- `v_id::Int` : node identifier
- `climate::Climate` : Climate data
- `calib_data::Array` : calibration data for target node by its id
- `metric::Function` : Optimization function to use. Defaults to RMSE.
"""
function calibrate!(sn::StreamfallNetwork, v_id::Int64, climate::Climate, calib_data::Array;
                    metric::Function=Streamfall.RMSE,
                    kwargs...)
    # Set defaults as necessary
    defaults = (;
        MaxTime=900,
        TraceInterval=30
    )
    kwargs = merge(defaults, kwargs)

    # Fitness of model is dependent on upstream node.
    ins = inlets(sn, v_id)

    # Recurse through and calibrate all nodes upstream
    if !isempty(ins)
        for nid in ins
            calibrate!(sn, nid, climate, calib_data; metric=metric, kwargs...)
        end
    end

    this_node = sn[v_id]

    # Create context-specific optimization function
    opt_func = x -> obj_func(x, climate, sn, v_id, calib_data; metric)

    # Get node parameters (default values and bounds)
    param_names, x0, param_bounds = param_info(this_node; with_level=false)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  kwargs...)

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated $(v_id) ($(this_node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(this_node, bs...)

    return res, opt
end


"""
    calibrate!(node, climate::Climate, calib_data::Union{Dict, DataFrame};
               metric::Function=Streamfall.RMSE,
               kwargs...)

Calibrate a given node using the BlackBoxOptim package.

# Arguments
- `node::NetworkNode` : Streamfall node
- `climate` : Climate data
- `calib_data` : calibration data for target node by its id
- `extractor::Function` : Calibration extraction method, define a custom one to change behavior
- `metric::Function` : Optimization function to use. Defaults to RMSE.
"""
function calibrate!(node::NetworkNode, climate::Climate, calib_data::Union{Dict, DataFrame};
                    metric::Function=Streamfall.RMSE,
                    kwargs...)
    # Set defaults as necessary
    defaults = (;
        MaxTime=900,
        TraceInterval=30
    )
    kwargs = merge(defaults, kwargs)

    # Create context-specific optimization function
    opt_func = x -> obj_func(x, climate, node, calib_data; metric)

    # Get node parameters (default values and bounds)
    param_names, x0, param_bounds = param_info(node; with_level=false)
    opt = bbsetup(opt_func; SearchRange=param_bounds,
                  kwargs...)

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated ($(node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(node, bs...)

    return res, opt
end


"""
    calibrate!(sn::StreamfallNetwork, climate::Climate, calib_data;
               metric::Function=Streamfall.RMSE,
               kwargs...)

Calibrate a stream network.
"""
function calibrate!(sn::StreamfallNetwork, climate::Climate, calib_data;
                   metric::Function=Streamfall.RMSE,
                   kwargs...)
    _, outlets = find_inlets_and_outlets(sn)
    calib_states = Dict()
    for out in outlets
        calib_states[out] = calibrate!(sn, out, climate, calib_data; metric=metric, kwargs...)
    end

    return calib_states
end


"""
Serialize calibration results and optimization object to disk.
"""
function save_calibration!(res, optobj, fn=nothing)
    if isnothing(fn)
        fn = "./temp" * string(rand(1:Int(1e8))) * ".tmp"
    end

    fh = open(fn, "w")
    serialize(fh, (res, optobj))
    close(fh)

    return fn
end


"""
Deserialize calibration results and optimization object from disk.
"""
function load_calibration(fn)
    fh = open(fn, "r")
    (res, optobj) = deserialize(fh)
    close(fh)

    return (res, optobj)
end
