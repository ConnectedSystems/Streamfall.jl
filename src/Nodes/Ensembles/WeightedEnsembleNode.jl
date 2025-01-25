Base.@kwdef mutable struct WeightedEnsemble{N<:NetworkNode, P, A<:Real} <: EnsembleNode
    name::String
    area::A

    instances::Array{N} = NetworkNode[]
    weights::Array{P} = [
        Param(1.0; bounds=[0.0,1.0])
    ]

    # Default to weighted sum
    comb_method::Function = (X, weights) -> sum([x * w for (x, w) in zip(X, weights)])

    outflow::Array{A} = []

    obj_func::Function = obj_func
end

function prep_state!(node::WeightedEnsemble, timesteps::Int64)::Nothing
    node.outflow = fill(0.0, timesteps)

    for n in node.instances
        prep_state!(n, timesteps)
    end

    return nothing
end


function param_info(node::WeightedEnsemble; kwargs...)::Tuple
    values = Float64[w.val for w in node.weights]
    bounds = [w.bounds for w in node.weights]
    param_names = ["weights_n$(i)" for i in 1:length(bounds)]

    return param_names, values, bounds
end


function WeightedEnsemble(nodes::Array{<:NetworkNode}; weights::Array{Float64},
                      bounds::Union{Nothing, Array{Tuple}}=nothing,
                      comb_method::Union{Nothing, Function}=nothing)::WeightedEnsemble
    if isnothing(bounds)
        num_nodes = length(nodes)
        @assert(length(weights) == num_nodes, "Number of nodes do not match provided number of weights")

        bounds = repeat([(0.0, 1.0)], num_nodes, 1)
    else
        @assert(length(bounds) == length(weights), "Number of bounds do not match provided number of weights")
    end

    p_weights = [Param(w; bounds=b) for (w, b) in zip(weights, bounds)]

    n1 = nodes[1]
    tmp = WeightedEnsemble{NetworkNode, Param, Float64}(;
        name=n1.name,
        area=n1.area,
        instances=nodes,
        weights=p_weights
    )

    if !isnothing(comb_method)
        tmp.comb_method = comb_method
    end

    return tmp
end


function run_node!(ensemble::WeightedEnsemble, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    for inst in ensemble.instances
        run_node!(inst, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    ensemble.outflow = ensemble.comb_method([inst.outflow for inst in ensemble.instances], ensemble.weights)
end


function run_timestep!(node::WeightedEnsemble, rain, et, ts; inflow=0.0, extraction=0.0, exchange=0.0)
    for inst in node.instances
        run_timestep!(inst, rain, et, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    node.outflow[ts] = node.comb_method([inst.outflow[ts] for inst in node.instances], node.weights)
end


function reset!(ensemble::WeightedEnsemble)
    for inst in ensemble.instances
        reset!(inst)
    end

    ensemble.outflow = []
end


# Functions to support calibration of weights
function update_params!(node::WeightedEnsemble, weights...)
    for (i, w) in enumerate(weights)
        node.weights[i] = Param(w, bounds=node.weights[i].bounds)
    end
end

"""
    calibrate_instances!(ensemble::WeightedEnsemble, climate::Climate, calib_data::Union{AbstractArray,DataFrame}, metric::Union{F,AbstractDict{String,F}}; kwargs...) where {F}

Calibrate individual model instances and the ensemble weights.

# Arguments
- `ensemble`: WeightedEnsemble node containing multiple model instances
- `climate`: Climate data for simulation
- `calib_data`: Calibration data, either as an array or DataFrame with node names as columns
- `metric`: Optimization metric function or Dict mapping node names to metrics
- `kwargs`: Additional arguments passed to BlackBoxOptim

# Returns
- Tuple of (optimization_result, optimization_setup) from weights calibration
"""
function calibrate_instances!(
    ensemble::WeightedEnsemble,
    climate::Climate,
    calib_data::Union{AbstractArray,DataFrame},
    metric::Union{F,AbstractDict{String,F}};
    kwargs...
) where {F}
    # Calibrate individual instances first
    for node in ensemble.instances
        calibrate!(node, climate, calib_data, metric; kwargs...)
    end

    # Then optimize weights
    return calibrate!(ensemble, climate, calib_data, metric; kwargs...)
end
