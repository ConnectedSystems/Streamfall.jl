using Parameters
using ModelParameters
using Statistics


abstract type EnsembleNode <: NetworkNode end

include("WeightedEnsembleNode.jl")
# include("StateEnsembleNode.jl")

