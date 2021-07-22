# Ensure dependent data and packages are available
using Distributed, BlackBoxOptim, Serialization


function data_extraction(node, calib_data::Dict)
    n_data = node.outflow
    h_data = calib_data[node.name]

    if h_data isa DataFrame
        h_data = Array(subcatchment_data(node, h_data, "_outflow"))
    end

    return h_data, n_data
end


function data_extraction(node, calib_data::DataFrame)
    n_data = node.outflow
    h_data = Array(subcatchment_data(node, calib_data, "_outflow"))

    return h_data, n_data
end


"""Calibrate current node."""
function obj_func(params, climate, sn, v_id, calib_data; extractor::Function, metric::Function)

    node = sn[v_id]
    update_params!(node, params...)

    ext = nothing
    fluxes = nothing
    try
        ext = calib_data[:, "$(node.name)_extractions"]
        fluxes = calib_data[:, "$(node.name)_exchange"]
    catch e
        if !isa(e, ArgumentError) && !isa(e, KeyError)
            throw(e)
        end
    end

    func! = get_prop(sn, v_id, :nfunc)
    func!(sn, v_id, climate; extraction=ext, exchange=fluxes)
    h_data, n_data = extractor(node, calib_data)
    score = metric(h_data, n_data)

    # Reset to clear stored values
    reset!(sn)

    return score
end


"""
    calibrate!(sn, v_id, climate, calib_data,
               extractor::Function=Streamfall.data_extraction,
               metric::Function=Streamfall.RMSE;
               kwargs...)

Calibrate a given node using the BlackBoxOptim package.

# Arguments
- `sn::StreamfallNetwork` : Network
- `v_id::Int` : node identifier
- `climate::Climate` : Climate data
- `calib_data::Union{Dict, DataFrame}` : calibration data for target node by its id
- `extractor::Function` : Calibration extraction method, define a custom one to change behavior
- `metric::Function` : Optimization function to use. Defaults to RMSE.
"""
function calibrate!(sn, v_id, climate, calib_data,
                   extractor::Function=Streamfall.data_extraction,
                   metric::Function=Streamfall.RMSE;
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
            calibrate!(sn, nid, climate, calib_data, extractor, metric; kwargs...)
        end
    end

    this_node = sn[v_id]

    # Create context-specific optimization function
    opt_func = x -> obj_func(x, climate, sn, v_id, calib_data; extractor, metric)

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
    calibrate!(sn, climate, calib_data,
               extractor::Function=Streamfall.data_extraction,
               metric::Function=Streamfall.RMSE;
               kwargs...)

Calibrate a stream network.
"""
function calibrate!(sn::StreamfallNetwork, climate::Climate, calib_data;
                   extractor::Function=Streamfall.data_extraction,
                   metric::Function=Streamfall.RMSE,
                   kwargs...)
    _, outlets = find_inlets_and_outlets(sn)
    calib_states = Dict()
    for out in outlets
        calib_states[out] = calibrate!(sn, out, climate, calib_data, extractor, metric; kwargs...)
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
