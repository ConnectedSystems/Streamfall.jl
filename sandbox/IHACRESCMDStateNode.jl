# IHACRES CMD State-based Node
using Streamfall
using ModelParameters, OnlineStats, Setfield
import Streamfall: run_node!, update_params!


Base.@kwdef mutable struct CMDStateNode{P, A<:Real} <: IHACRESNode
    Streamfall.@network_node

    # Activate parameters
    d::A = 200.0;
    d2::A = 2.0
    e::A = 1.0
    f::A = 0.8
    a::A = 0.9
    b::A = 0.1
    storage_coef::A = 2.9
    alpha::A = 0.1

    # low CMD (wet)
    low_d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    low_d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    low_e::P = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    low_f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    low_a::P = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == (1/tau_q)
    low_b::P = Param(0.1, bounds=(1e-3, 0.1))  # slowflow storage coefficent == (1/tau_s)
    low_storage_coef::P = Param(2.9, bounds=(1e-10, 10.0))
    low_alpha::P = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    # mid CMD ("usual")
    mid_d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    mid_d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    mid_e::P = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    mid_f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    mid_a::P = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == (1/tau_q)
    mid_b::P = Param(0.1, bounds=(1e-3, 0.1))  # slowflow storage coefficent == (1/tau_s)
    mid_storage_coef::P = Param(2.9, bounds=(1e-10, 10.0))
    mid_alpha::P = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    # high CMD ("dry")
    high_d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    high_d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    high_e::P = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    high_f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    high_a::P = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == (1/tau_q)
    high_b::P = Param(0.1, bounds=(1e-3, 0.1))  # slowflow storage coefficent == (1/tau_s)
    high_storage_coef::P = Param(2.9, bounds=(1e-10, 10.0))
    high_alpha::P = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    level_params::Array{P, 1} = [
        Param(-0.01, bounds=(-10.0, -0.01)),  # p1
        Param(0.8, bounds=(0.0, 1.5)),  # p2
        Param(4.5, bounds=(0.0, 20.0)), # p3
        Param(5.0, bounds=(1.0, 10.0)), # p4
        Param(0.35, bounds=(0.0, 1.0)), # p5
        Param(1.41, bounds=(-2.0, 2.0)), # p6
        Param(-1.45, bounds=(-2.5, 0.0)), # p7
        Param(6.75, bounds=(0.0, 10.0)), # p8
        Param(150.0, bounds=(50.0, 200.0)) # ctf
    ]

    storage::Array{A} = [100.0]  # CMD
    quick_store::Array{A} = [0.0]
    slow_store::Array{A} = [0.0]
    outflow::Array{A} = []
    effective_rainfall::Array{A} = []
    et::Array{A} = []
    inflow::Array{A} = []
    level::Array{A} = []
    gw_store::Array{A} = [0.0]

    # cache arrays
    cache::NamedTuple{(:i_cache, :r_cache), <:Tuple{Vector{Float64}, Vector{Float64}}} = (i_cache=zeros(3), r_cache=zeros(2))

    # active param records
    active_param_set::Array{Int64} = [1]
end


"""Generic run_node method that handles any single state variable.

Develops an "online" running statistic of states from the first N days, 
which remain static for the rest of the analysis period.
"""
function run_node!(node::CMDStateNode, climate::Climate; quantiles::Array=[0.1, 0.9], burn_in=1826, releases=nothing, log_transform=false)::Nothing

    # quantiles = [0.0, 0.1, 0.9]
    # new quantiles = [0.1, 0.9]

    # Get node parameters
    _, x0, __ = param_info(node; with_level=false)
    # num_params = length(x0)

    timesteps = length(climate)
    active_param_set = zeros(Int, timesteps)

    # state_var = getfield(node, state)
    state_var = node.storage

    o = Quantile(quantiles, b=2000)
    param_idxs = nothing
    thresholds = nothing
    for ts in (1:timesteps)
        # Update param set based on state
        state_val = state_var[ts]
        if log_transform
            # add constant to avoid log(0)
            tmp = log(state_val .+ 1e-2)
        else
            tmp = state_val
        end

        if ts < burn_in
            # Update thresholds based on state value so far
            fit!(o, tmp)
            thresholds = value(o)
        end

        active_param_set[ts] = switch_params!(node, tmp, thresholds)

        # param_set_id, node_params, param_idxs = find_state_vars(tmp, thresholds, params, num_params, length(thresholds))
        # if (ts == 1) || (param_set_id != active_param_set[ts-1])
        #     update_params!(node, node_params...)
        # end

        # record timestep in which this param was active
        # active_param_set[ts] = param_set_id

        Streamfall.run_node!(node, climate, ts; extraction=releases)
    end

    node.active_param_set = active_param_set

    # return active_param_set, param_idxs
    return nothing
end


"""
Clunky approach to switching between parameter sets
"""
function switch_params!(node::CMDStateNode, state_val::T, thresholds::Vector{T})::Int64 where {T<:Real}
    set_id = first(searchsortedfirst.(Ref(thresholds), state_val))
    if set_id == node.active_param_set[end]
        return set_id
    end

    if set_id == 1
        node.d = node.low_d.val
        node.d2 = node.low_d2.val
        node.e = node.low_e.val
        node.f = node.low_f.val
        node.a = node.low_a.val
        node.b = node.low_b.val
        node.storage_coef = node.low_storage_coef.val
        node.alpha = node.low_alpha.val
    elseif set_id == 2
        node.d = node.mid_d.val
        node.d2 = node.mid_d2.val
        node.e = node.mid_e.val
        node.f = node.mid_f.val
        node.a = node.mid_a.val
        node.b = node.mid_b.val
        node.storage_coef = node.mid_storage_coef.val
        node.alpha = node.mid_alpha.val
    elseif set_id == 3
        node.d = node.high_d.val
        node.d2 = node.high_d2.val
        node.e = node.high_e.val
        node.f = node.high_f.val
        node.a = node.high_a.val
        node.b = node.high_b.val
        node.storage_coef = node.high_storage_coef.val
        node.alpha = node.high_alpha.val
    else
        error("Unknown state id")
    end

    return set_id
end


"""

Clunky implementation of updating parameters.
On update, 
"""
function update_params!(node::CMDStateNode, params...)::Nothing
    node.d = params[1]
    node.d2 = params[2]
    node.e = params[3]
    node.f = params[4]
    node.a = params[5]
    node.b = params[6]
    node.storage_coef = params[7]
    node.alpha = params[8]

    node.low_d = Param(params[1], bounds=node.low_d.bounds::Tuple)
    node.low_d2 = Param(params[2], bounds=node.low_d2.bounds::Tuple)
    node.low_e = Param(params[3], bounds=node.low_e.bounds::Tuple)
    node.low_f = Param(params[4], bounds=node.low_f.bounds::Tuple)
    node.low_a = Param(params[5], bounds=node.low_a.bounds::Tuple)
    node.low_b = Param(params[6], bounds=node.low_b.bounds::Tuple)
    node.low_storage_coef = Param(params[7], bounds=node.low_storage_coef.bounds::Tuple)
    node.low_alpha = Param(params[8], bounds=node.low_alpha.bounds::Tuple)

    node.mid_d = Param(params[9], bounds=node.mid_d.bounds::Tuple)
    node.mid_d2 = Param(params[10], bounds=node.mid_d2.bounds::Tuple)
    node.mid_e = Param(params[11], bounds=node.mid_e.bounds::Tuple)
    node.mid_f = Param(params[12], bounds=node.mid_f.bounds::Tuple)
    node.mid_a = Param(params[13], bounds=node.mid_a.bounds::Tuple)
    node.mid_b = Param(params[14], bounds=node.mid_b.bounds::Tuple)
    node.mid_storage_coef = Param(params[15], bounds=node.mid_storage_coef.bounds::Tuple)
    node.mid_alpha = Param(params[16], bounds=node.mid_alpha.bounds::Tuple)

    node.high_d = Param(params[17], bounds=node.high_d.bounds::Tuple)
    node.high_d2 = Param(params[18], bounds=node.high_d2.bounds::Tuple)
    node.high_e = Param(params[19], bounds=node.high_e.bounds::Tuple)
    node.high_f = Param(params[20], bounds=node.high_f.bounds::Tuple)
    node.high_a = Param(params[21], bounds=node.high_a.bounds::Tuple)
    node.high_b = Param(params[22], bounds=node.high_b.bounds::Tuple)
    node.high_storage_coef = Param(params[23], bounds=node.high_storage_coef.bounds::Tuple)
    node.high_alpha = Param(params[24], bounds=node.high_alpha.bounds::Tuple)

    # Below doesn't seem to work yet...
    # @set node.low_d.val = params[1]
    # @set node.low_d2.val = params[2]
    # @set node.low_e.val = params[3]
    # @set node.low_f.val = params[4]
    # @set node.low_a.val = params[5]
    # @set node.low_b.val = params[6]
    # @set node.low_storage_coef.val = params[7]
    # @set node.low_alpha.val = params[8]

    # @set node.mid_d.val = params[9]
    # @set node.mid_d2.val = params[10]
    # @set node.mid_e.val = params[11]
    # @set node.mid_f.val = params[12]
    # @set node.mid_a.val = params[13]
    # @set node.mid_b.val = params[14]
    # @set node.mid_storage_coef.val = params[15]
    # @set node.mid_alpha.val = params[16]

    # @set node.high_d.val = params[17]
    # @set node.high_d2.val = params[18]
    # @set node.high_e.val = params[19]
    # @set node.high_f.val = params[20]
    # @set node.high_a.val = params[21]
    # @set node.high_b.val = params[22]
    # @set node.high_storage_coef.val = params[23]
    # @set node.high_alpha.val = params[24]

    node.active_param_set = [1]

    return nothing
end


function run_node!(node::CMDStateNode, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    invoke(run_node!, Tuple{IHACRESNode, Climate}, node, climate; inflow, extraction, exchange)
end


# function find_state_vars(val, thresholds, params, num_params, n_states)
#     # index of parameters (e.g., if 3 partitions, then create 3 arrays of 8 parameters)
#     param_idxs = collect(Iterators.partition(1:length(params), num_params))

#     # searchsortedfirst.(Ref([3, 5, 10]), [2, 2, 2, 12])

#     # TODO: Check if we need to exclude first bin
#     set_id = searchsortedfirst.(Ref(thresholds), val)



#     param_set = param_idxs[set_id]
#     node_params = params[param_set]

#     return set_id, node_params, param_idxs
# end
