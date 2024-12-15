using Serialization
using BlackBoxOptim
using Distributed


"""Calibrate the specified node in network."""
function obj_func(
    params, climate::Climate, sn::StreamfallNetwork, v_id::Int, calib_data::Array;
    metric::Function, inflow=nothing, extraction=nothing, exchange=nothing
)
    return obj_func(params, climate, sn[v_id], calib_data; metric=metric, inflow=inflow, extraction=extraction, exchange=exchange)
end


"""Calibrate current node."""
function obj_func(
    params, climate::Climate, node::NetworkNode, calib_data::Array;
    metric::Function, inflow=nothing, extraction=nothing, exchange=nothing
)
    update_params!(node, params...)

    run_node!(node, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    score = metric(calib_data, node.outflow)

    # Reset to clear stored values
    reset!(node)

    return score
end


"""
    dependent_obj_func(
        params, climate::Climate, this_node::NetworkNode, next_node::NetworkNode, calib_data::DataFrame;
        metric::Function, weighting=0.5, inflow=nothing, extraction=nothing, exchange=nothing
    )

Objective function which considers performance of next node and current node.
The weighting factor places equal emphasis on both nodes by default (0.5).
The weighting value places a \$x\$ weight on the current node, and \$1 - x\$ on the next
node.

"""
function dependent_obj_func(
    params, climate::Climate, this_node::NetworkNode, next_node::NetworkNode, calib_data::DataFrame;
    metric::Function, weighting=0.5, inflow=nothing, extraction=nothing, exchange=nothing
)
    update_params!(this_node, params...)

    # Run dependent nodes
    run_node!(
        this_node, climate;
        inflow=inflow, extraction=extraction, exchange=exchange
    )
    run_node!(
        next_node, climate;
        inflow=this_node.outflow, extraction=extraction, exchange=exchange
    )

    # Alias data as necessary
    if typeof(next_node) <: DamNode
        # Calibrate against both outflow and dam levels
        sim_data = next_node.level
        obs_data = calib_data[:, next_node.name]

        if weighting == 0.0
            # All emphasis is on dam levels
            score = metric(obs_data, sim_data)
        elseif weighting < 1.0
            # Use a mix
            dam_score = metric(obs_data, sim_data)

            sim_data = this_node.outflow
            obs_data = calib_data[:, this_node.name]

            flow_score = metric(obs_data, sim_data)

            # Use weighted average of the two
            score = (flow_score * weighting) + (dam_score * (1.0 - weighting))
        else
            # Use just outflows
            sim_data = this_node.outflow
            obs_data = calib_data[:, this_node.name]

            score = metric(obs_data, sim_data)
        end
    elseif typeof(this_node) <: DamNode
        # Calibrate against dam levels only
        sim_data = this_node.level
        obs_data = calib_data[:, this_node.name]

        score = metric(obs_data, sim_data)
    else
        # Calibrate against outflows only
        sim_data = this_node.outflow
        obs_data = calib_data[:, this_node.name]

        score = metric(obs_data, sim_data)
    end

    reset!(this_node)
    reset!(next_node)

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
    calibrate!(sn, v_id, climate, calib_data; metric=metric, kwargs...)
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
- `calib_data::DataFrame` : calibration data for target node by its id
- `metric::Function` : Optimization function to use. Defaults to RMSE.
"""
function calibrate!(sn::StreamfallNetwork, v_id::Int64, climate::Climate, calib_data::DataFrame;
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
    next_node = sn[first(outlets(sn, this_node.name))]

    extraction = get(kwargs, :extraction, nothing)
    exchange = get(kwargs, :exchange, nothing)

    # Create context-specific optimization function
    if typeof(next_node) <: DamNode
        # If the next node represents a dam, attempt to calibrate considering the outflows
        # from the current node and the Dam Levels of the next node.
        weighting = get(kwargs, :weighting, 0.5)
        opt_func = x -> next_node.obj_func(
            x, climate, this_node, next_node, calib_data;
            metric=metric, extraction=extraction, exchange=exchange, weighting=weighting
        )
    elseif typeof(this_node) <: DamNode
        opt_func = x -> obj_func(
            x, climate, this_node, calib_data[:, this_node.name];
            metric=metric, extraction=extraction, exchange=exchange
        )
    else
        opt_func = x -> this_node.obj_func(
            x, climate, this_node, calib_data[:, this_node.name];
            metric=metric, extraction=extraction, exchange=exchange
        )
    end

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


# TODO : Clean the next two methods up as they are rough duplicates.
"""
    calibrate!(node, climate::Climate, calib_data::Array;
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
function calibrate!(node::NetworkNode, climate::Climate, calib_data::Array;
                    metric::Function=Streamfall.RMSE,
                    kwargs...)
    # Set defaults as necessary
    defaults = (;
        MaxTime=900,
        TraceInterval=30
    )
    kwargs = merge(defaults, kwargs)

    next_node = sn[first(outlets(sn, node.name))]

    extraction = get(kwargs, :extraction, nothing)
    exchange = get(kwargs, :exchange, nothing)

    # Create context-specific optimization function
    if next_node <: DamNode
        # If the next node represents a dam, attempt to calibrate considering the outflows
        # from the current node and the Dam Levels of the next node.
        opt_func = x -> next_node.obj_func(x, climate, node, next_node, calib_data; metric=metric, extraction=extraction, exchange=exchange)
    else
        opt_func = x -> this_node.obj_func(x, climate, node, calib_data; metric=metric, extraction=extraction, exchange=exchange)
    end

    # Get node parameters (default values and bounds)
    param_names, x0, param_bounds = param_info(node; with_level=false)
    opt = bbsetup(
        opt_func;
        SearchRange=param_bounds,
        x0=x0,
        kwargs...
    )

    res = bboptimize(opt)

    bs = best_candidate(res)
    @info "Calibrated ($(node.name)), with score: $(best_fitness(res))"
    @info "Best Params:" collect(bs)

    # Update node with calibrated parameters
    update_params!(node, bs...)

    return res, opt
end



"""
    calibrate!(node, climate::Climate, calib_data::DataFrame;
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
function calibrate!(node::NetworkNode, climate::Climate, calib_data::DataFrame;
                    metric::Function=Streamfall.RMSE,
                    kwargs...)
    # Set defaults as necessary
    defaults = (;
        MaxTime=900,
        TraceInterval=30
    )
    kwargs = merge(defaults, kwargs)

    # reset!(node)

    # Create context-specific optimization function
    opt_func = x -> obj_func(x, climate, node, calib_data[:, node.name]; metric=metric)

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
    calibrate!(sn::StreamfallNetwork, climate::Climate, calib_data::DataFrame;
               metric::Function=Streamfall.RMSE,
               kwargs...)

Calibrate a stream network.
"""
function calibrate!(sn::StreamfallNetwork, climate::Climate, calib_data::DataFrame;
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
