using Parameters
using ModelParameters
using Statistics


abstract type EnsembleNode <: NetworkNode end

include("WeightedEnsembleNode.jl")
include("GREnsembleNode.jl")
# include("StateEnsembleNode.jl")

function prep_state!(node::EnsembleNode, timesteps::Int64)::Nothing
    node.outflow = fill(0.0, timesteps)

    for n in node.instances
        prep_state!(n, timesteps)
    end

    return nothing
end
