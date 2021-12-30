using Parameters
using ModelParameters
using Statistics


abstract type EnsembleNode <: NetworkNode end


Base.@kwdef mutable struct BaseEnsemble{N<:NetworkNode, P, R<:Real} <: EnsembleNode
    @network_node

    instances::Array{N} = []
    weights::Array{P} = [
        Param(1.0; bounds=[0.0,1.0])
    ]

    # Default to weighted sum
    comb_method::Function = (X, weights) -> sum([x * w for (x, w) in zip(X, weights)])

    outflow::Array{R} = []
end


function param_info(node::BaseEnsemble; kwargs...)::Tuple
    values = Float64[w.val for w in node.weights]
    bounds = [w.bounds for w in node.weights]
    param_names = ["weights_n$(i)" for i in 1:length(bounds)]

    return param_names, values, bounds
end


function BaseEnsemble(nodes::Array{<:NetworkNode}; weights::Array{Float64},
                      bounds::Union{Nothing, Array{Tuple}}=nothing,
                      comb_method::Union{Nothing, Function}=nothing)::BaseEnsemble
    if isnothing(bounds)
        num_nodes = length(nodes)
        @assert(length(weights) == num_nodes, "Number of nodes do not match provided number of weights")

        bounds = repeat([(0.0, 1.0)], num_nodes, 1)
    else
        @assert(length(bounds) == length(weights), "Number of bounds do not match provided number of weights")
    end

    p_weights = [Param(w; bounds=b) for (w, b) in zip(weights, bounds)]

    n1 = nodes[1]
    if !isnothing(comb_method)
        tmp = BaseEnsemble{NetworkNode, Param, Float64}(; name=n1.name, area=n1.area, instances=nodes, weights=p_weights, comb_method)
    else
        tmp = BaseEnsemble{NetworkNode, Param, Float64}(; name=n1.name, area=n1.area, instances=nodes, weights=p_weights)
    end

    return tmp
end


function run_node!(ensemble::BaseEnsemble, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    for inst in ensemble.instances
        run_node!(inst, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    ensemble.outflow = ensemble.comb_method([inst.outflow for inst in ensemble.instances], ensemble.weights)
end


function run_timestep!(node, rain, et, ts; inflow=0.0, extraction=0.0, exchange=0.0)
    for inst in node.instances
        run_timestep!(inst, rain, et, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    Qt = node.comb_method([inst.outflow[ts] for inst in node.instances], node.weights)

    append!(node.outflow, Qt)
end


function reset!(ensemble::BaseEnsemble)
    for inst in ensemble.instances
        reset!(inst)
    end

    ensemble.outflow = []
end


# Functions to support calibration of weights
function update_params!(node::BaseEnsemble, weights...)
    for (i, w) in enumerate(weights)
        node.weights[i] = Param(w, bounds=node.weights[i].bounds)
    end
end


# Calibrate each individual instance?
# function calibrate!(ensemble::BaseEnsemble, climate::Climate)
#     ensemble.weights
# end
