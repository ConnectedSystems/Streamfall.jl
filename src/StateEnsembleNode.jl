# using Parameters
# using ModelParameters


# Base.@kwdef mutable struct StateEnsemble <: EnsembleNode
#     instances::Array
#     comb_method::Function

#     outflow::Array = []
# end


# function BaseEnsemble(nodes::T...) where T <: NetworkNode
#     return BaseEnsemble([nodes...])
# end


# function run_node!(ensemble::StateEnsemble, climate::Climate)
#     for inst in ensemble.instance
#         run_node!(inst, climate)
#     end

#     ensemble.outflow = comb_method([inst.outflow for inst in ensemble.instance]...)
# end


# function reset!(ensemble::StateEnsemble)
#     for inst in ensemble.instance
#         reset!(inst)
#     end

#     ensemble.outflow = []
# end