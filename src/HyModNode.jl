using Parameters
using ModelParameters


abstract type HyModNode <: NetworkNode end


"""
Simple implementation of HyMod - does not include snow melt processes (see [1]).

Adapted with kind permission from: https://github.com/jdherman/GRA-2020-SALib

# References
1. Gharari, S., Hrachowitz, M., Fenicia, F., Savenije, H.H.G., 2013. 
    An approach to identify time consistent model parameters: sub-period calibration. 
    Hydrology and Earth System Sciences 17, 149–161. 
    https://doi.org/10.5194/hess-17-149-2013
"""
Base.@kwdef mutable struct SimpleHyModNode{P} <: HyModNode
    @network_node

    # parameters
    Sm_max::P = Param(250.0, bounds=(1.0, 500.0))
    B::P = Param(1.0, bounds=(0.0, 2.0))
    alpha::P = Param(0.2, bounds=(0.0, 1.0))
    Kf::P = Param(0.5, bounds=(0.1, 0.9999))
    Ks::P = Param(0.05, bounds=(0.001, 0.1))

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


# function prep_state!(node::HyModNode, sim_length::Int64)
#     node.Sm = fill(0.0, sim_length+1)
#     node.Sf1 = fill(0.0, sim_length+1)
#     node.Sf3 = fill(0.0, sim_length+1)
#     node.Ss1 = fill(0.0, sim_length+1)

#     node.outflow = fill(0.0, sim_length)
# end


"""
    run_node!(node::HyModNode, climate::Climate;
              inflow=nothing, extraction=nothing, exchange=nothing)

Run given HyMod node for entire simulation period.
"""
function run_node!(node::HyModNode, climate::Climate;
                   inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    # prep_state!(node, timesteps)

    for ts in 1:timesteps
        run_timestep!(node, climate, ts;
                      inflow=inflow, extraction=extraction, exchange=exchange)
    end
end


"""
    run_timestep!(node::HyModNode, climate::Climate, timestep::Int,
                  inflow::Float64, extraction::Float64, exchange::Float64)

Run given HyMod node for a time step.
"""
function run_timestep!(node::SimpleHyModNode, climate::Climate, timestep::Int;
                       inflow=nothing, extraction=nothing, exchange=nothing)::Float64
    ts = timestep
    P, ET = climate_values(node, climate, ts)
    return run_timestep!(node, P, ET, ts; inflow=inflow, extraction=extraction, exchange=exchange)
end


function run_timestep!(node::SimpleHyModNode, rain, et, ts; inflow=nothing, extraction=nothing, exchange=nothing)
    in_flow = timestep_value(ts, node.name, "inflow", inflow)
    ext = timestep_value(ts, node.name, "extraction", extraction)
    flux = timestep_value(ts, node.name, "exchange", exchange)

    Sm = node.Sm[ts]
    Sf1 = node.Sf1[ts]
    Sf2 = node.Sf2[ts]
    Sf3 = node.Sf3[ts]
    Ss1 = node.Ss1[ts]

    return run_hymod!(node, rain, et, Sm, Sf1, Sf2, Sf3, Ss1, in_flow, ext, flux)
end


function run_hymod!(node::SimpleHyModNode, P, PET, Sm, Sf1, Sf2, Sf3, Ss1, inflow, ext, flux)
    Sm_max = node.Sm_max
    B = node.B
    alpha = node.alpha
    Kf = node.Kf
    Ks = node.Ks

    # Calculate fluxes
    Peff = P*(1 - max(1.0 - Sm/Sm_max, 0.0)^B)  # PDM model Moore 1985
    ETa = min(PET*(Sm/Sm_max), Sm)

    Qf1 = Kf*Sf1
    Qf2 = Kf*Sf2
    Qf3 = Kf*Sf3
    Qs1 = Ks*Ss1

    # update state variables
    Sm_t1 = max(Sm + P - Peff - ETa, 0.0)
    Sf1_t1 = Sf1 + alpha*Peff - Qf1
    Sf2_t1 = Sf2 + Qf1 - Qf2
    Sf3_t1 = Sf3 + Qf2 - Qf3
    Ss1_t1 = Ss1 + (1-alpha)*Peff - Qs1

    Qtmp = (Qs1 + Qf3) * node.area
    Q_t1 = inflow + Qtmp - ext + flux
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
