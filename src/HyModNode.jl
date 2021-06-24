using Parameters
using ModelParameters


abstract type HyModNode <: NetworkNode end


"""
Simple implementation of HyMod - does not include snow melt processes.
"""
Base.@kwdef mutable struct SimpleHyModNode{A <: Union{Param, Real}} <: HyModNode
    @network_node

    # parameters
    Sm_max::A = Param(250.0, bounds=(1.0, 500.0))
    B::A = Param(1.0, bounds=(0.1, 2.0))
    alpha::A = Param(0.2, bounds=(0.1, 0.999))
    Kf::A = Param(0.5, bounds=(0.1, 0.999))
    Ks::A = Param(0.05, bounds=(0.001, 0.1))

    # stores
    Sm::Array{Float64} = [0.0]
    Sf1::Array{Float64} = [0.0]
    Sf2::Array{Float64} = [0.0]
    Sf3::Array{Float64} = [0.0]
    Ss1::Array{Float64} = [0.0]

    outflow::Array{Float64} = []
end


function SimpleHyModNode(name::String, spec::Dict)
    n = SimpleHyModNode{Param}(; name=name, area=spec["area"])
    node_params = spec["parameters"]
    for (p_name, p_val) in node_params
        sym = Symbol(p_name)
        p = getfield(n, sym)
        setfield!(n, sym, Param(p_val, bounds=p.bounds))
    end

    return n
end


function SimpleHyModNode(name::String, area::Float64, sm_max::Float64, B::Float64,
                         alpha::Float64, Kf::Float64, Ks::Float64)
    n = SimpleHyModNode{Param}(; name=name, area=area)
    update_params!(n, sm_max, B, alpha, Kf, Ks)

    # stores/output
    n.Sm = Float64[0.0]
    n.Sf1 = Float64[0.0]
    n.Sf2 = Float64[0.0]
    n.Sf3 = Float64[0.0]
    n.Ss1 = Float64[0.0]
    n.outflow = Float64[]

    return n
end


function run_node!(node::SimpleHyModNode, climate; extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    for ts in 1:timesteps
        run_node!(node, climate, ts; extraction=extraction, exchange=exchange)
    end

    return node.outflow
end


function run_node!(node::HyModNode, climate, ts; extraction=nothing, exchange=nothing)::Float64
    P, ET = climate_values(node, climate, ts)
    ext = timestep_value(ts, node.name, "extraction", extraction)
    flux = timestep_value(ts, node.name, "exchange", exchange)

    Sm = node.Sm[ts]
    Sm_max = node.Sm_max
    B = node.B
    alpha = node.alpha
    Kf = node.Kf
    Ks = node.Ks

    Sf1 = node.Sf1[ts]
    Sf2 = node.Sf2[ts]
    Sf3 = node.Sf3[ts]
    Ss1 = node.Ss1[ts]

    tmp_Q = run_hymod(node, P, ET, Sm, Sm_max, B, alpha, Kf, Ks, Sf1, Sf2, Sf3, Ss1)
    Q_t1 = tmp_Q - ext + flux

    update_states(node, Sm_t1, Sf1_t1, Sf2_t1, Sf3_t1, Ss1_t1, Q_t1)

    return Q_t1
end


"""
    run_node!(node::HyModNode, rain::Float64, evap::Float64, 
              inflow::Float64, extraction::Float64, exchange::Float64,
              timestep::Int)

Run given HyMod node for a time step (or last known state if not provided).
"""
function run_node!(node::HyModNode, rain::Float64, evap::Float64, 
                   inflow::Float64, extraction::Float64, exchange::Float64, timestep=nothing)
    ext = extraction
    flux = exchange
    ts = timestep
    if isnothing(ts)
        Sm = node.Sm[end]
        Sf1 = node.Sf1[end]
        Sf2 = node.Sf2[end]
        Sf3 = node.Sf3[end]
        Ss1 = node.Ss1[end]
    else
        Sm = node.Sm[ts]
        Sf1 = node.Sf1[ts]
        Sf2 = node.Sf2[ts]
        Sf3 = node.Sf3[ts]
        Ss1 = node.Ss1[ts]
    end

    Sm_max = node.Sm_max
    B = node.B
    alpha = node.alpha
    Kf = node.Kf
    Ks = node.Ks

    tmp_Q = run_hymod(node, rain, evap, Sm, Sm_max, B, alpha, Kf, Ks, Sf1, Sf2, Sf3, Ss1)
    Q_t1 = inflow + tmp_Q - ext + flux

    return Q_t1
end


function run_hymod(node::HyModNode, P, ET, Sm, Sm_max, B, alpha, Kf, Ks, Sf1, Sf2, Sf3, Ss1)
    # Calculate fluxes
    Peff = P*(1 - max(1.0 - Sm/Sm_max,0.0)^B) # PDM model Moore 1985
    Evap = min(ET*(Sm/Sm_max), Sm)

    Qf1 = Kf*Sf1
    Qf2 = Kf*Sf2
    Qf3 = Kf*Sf3
    Qs1 = Ks*Ss1

    # update state variables
    Sm_t1 = Sm + P - Peff - Evap
    Sf1_t1 = Sf1 + alpha*Peff - Qf1
    Sf2_t1 = Sf2 + Qf1 - Qf2
    Sf3_t1 = Sf3 + Qf2 - Qf3
    Ss1_t1 = Ss1 + (1-alpha)*Peff - Qs1

    Q_t1 = (Qs1 + Qf3)

    update_states(node, Sm_t1, Sf1_t1, Sf2_t1, Sf3_t1, Ss1_t1, Q_t1)

    return Q_t1
end


function update_params!(node::HyModNode, Sm_max, B, alpha, Kf, Ks)
    node.Sm_max = Param(Sm_max, bounds=node.Sm_max.bounds::Tuple)
    node.B = Param(B, bounds=node.B.bounds::Tuple)
    node.alpha = Param(alpha, bounds=node.alpha.bounds::Tuple)
    node.Kf = Param(Kf, bounds=node.Kf.bounds::Tuple)
    node.Ks = Param(Ks, bounds=node.Ks.bounds::Tuple)
end


function update_states(node::HyModNode, Sm, Sf1, Sf2, Sf3, Ss1, Q)
    append!(node.Sm, Sm)
    append!(node.Sf1, Sf1)
    append!(node.Sf2, Sf2)
    append!(node.Sf3, Sf3)
    append!(node.Ss1, Ss1)
    append!(node.outflow, Q)
end


function reset!(node::HyModNode)
    node.Sm = Float64[node.Sm[1]]
    node.Sf1 = Float64[node.Sf1[1]]
    node.Sf2 = Float64[node.Sf2[1]]
    node.Sf3 = Float64[node.Sf3[1]]
    node.Ss1 = Float64[node.Ss1[1]]
    node.outflow = Float64[]
end
