using Parameters
using ModelParameters


abstract type HyModNode <: NetworkNode end


"""
Simple implementation of HyMod - does not include snow melt processes (see [1]).

Adapted with kind permission from: https://github.com/jdherman/GRA-2020-SALib

# References
1. Gharari, S., Hrachowitz, M., Fenicia, F., Savenije, H.H.G., 2013.
    An approach to identify time consistent model parameters: sub-period calibration.
    Hydrology and Earth System Sciences 17, 149-161.
    https://doi.org/10.5194/hess-17-149-2013
2. Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, H. V.,
    and Sorooshian, S., 2001. A framework for development and applica-
    tion of hydrological models, Hydrol. Earth Syst. Sci., 5, 13-26,
    https://doi.org/10.5194/hess-5-13-2001.
"""
Base.@kwdef mutable struct SimpleHyModNode{P,A<:AbstractFloat} <: HyModNode
    const name::String
    const area::A

    # parameters
    Sm_max::P = Param(250.0, bounds=(1.0, 500.0))
    B::P = Param(1.0, bounds=(0.0, 2.0))
    alpha::P = Param(0.2, bounds=(0.0, 1.0))
    Kf::P = Param(0.5, bounds=(0.1, 0.9999))
    Ks::P = Param(0.05, bounds=(0.001, 0.1))

    # stores
    Sm::Array{A} = [0.0]
    Sf1::Array{A} = [0.0]
    Sf2::Array{A} = [0.0]
    Sf3::Array{A} = [0.0]
    Ss1::Array{A} = [0.0]

    outflow::Array{A} = []

    obj_func::Function = obj_func
end


function SimpleHyModNode(name::String, spec::AbstractDict)
    n = create_node(SimpleHyModNode, name, spec["area"])
    node_params = spec["parameters"]
    for (p_name, p_val) in node_params
        sym = Symbol(p_name)
        p = getfield(n, sym)
        setfield!(n, sym, Param(p_val, bounds=p.bounds))
    end

    return n
end

"""
    SimpleHyModNode(
        name::String, area::Float64, sm_max::Float64, B::Float64,
        alpha::Float64, Kf::Float64, Ks::Float64
    )
"""
function SimpleHyModNode(
    name::String, area::Float64, sm_max::Float64, B::Float64,
    alpha::Float64, Kf::Float64, Ks::Float64
)
    n = create_node(SimpleHyModNode, name, area)
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


function prep_state!(node::SimpleHyModNode, sim_length::Int64)::Nothing
    resize!(node.Sm, sim_length + 1)
    resize!(node.Sf1, sim_length + 1)
    resize!(node.Sf2, sim_length + 1)
    resize!(node.Sf3, sim_length + 1)
    resize!(node.Ss1, sim_length + 1)

    node.Sm[2:end] .= 0.0
    node.Sf1[2:end] .= 0.0
    node.Sf2[2:end] .= 0.0
    node.Sf3[2:end] .= 0.0
    node.Ss1[2:end] .= 0.0

    resize!(node.outflow, sim_length)
    node.outflow .= 0.0

    return nothing
end

"""
    run_timestep!(
        node::SimpleHyModNode, climate::Climate, timestep::Int;
        inflow=nothing, extraction=nothing, exchange=nothing
    )::Float64

Run given HyMod node for a time step.
"""
function run_timestep!(
    node::SimpleHyModNode, climate::Climate, timestep::Int;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Float64
    ts = timestep
    P, ET = climate_values(node, climate, ts)

    return run_timestep!(node, P, ET, ts; inflow=inflow, extraction=extraction, exchange=exchange)
end


function run_timestep!(
    node::SimpleHyModNode, rain::F, et::F, ts::Int64;
    inflow=nothing, extraction=nothing, exchange=nothing
)::F where {F<:Float64}
    in_flow = timestep_value(ts, node.name, "inflow", inflow)
    ext = timestep_value(ts, node.name, "extraction", extraction)
    flux = timestep_value(ts, node.name, "exchange", exchange)

    return run_hymod!(node, ts, rain, et, in_flow, ext, flux)
end


function run_hymod!(node::SimpleHyModNode, ts::Int64, P::F, PET::F, inflow::F, ext::F, flux::F) where {F<:Float64}
    Sm::F = node.Sm[ts]
    Sf1::F = node.Sf1[ts]
    Sf2::F = node.Sf2[ts]
    Sf3::F = node.Sf3[ts]
    Ss1::F = node.Ss1[ts]

    Sm_max::F = node.Sm_max.val::F
    B::F = node.B.val::F
    alpha::F = node.alpha.val::F
    Kf::F = node.Kf.val::F
    Ks::F = node.Ks.val::F

    # Calculate fluxes
    Peff::F = P * (1.0 - max(1.0 - Sm / Sm_max, 0.0)^B)  # PDM model Moore 1985
    ETa::F = min(PET * (Sm / Sm_max), Sm)

    Qf1::F = Kf * Sf1
    Qf2::F = Kf * Sf2
    Qf3::F = Kf * Sf3
    Qs1::F = Ks * Ss1

    # update state variables
    Sm_t1::F = max(Sm + P - Peff - ETa, 0.0)
    Sf1_t1::F = Sf1 + alpha * Peff - Qf1
    Sf2_t1::F = Sf2 + Qf1 - Qf2
    Sf3_t1::F = Sf3 + Qf2 - Qf3
    Ss1_t1::F = Ss1 + (1 - alpha) * Peff - Qs1

    Qtmp::F = (Qs1 + Qf3) * node.area
    Q_t::F = inflow + Qtmp - ext + flux
    update_state!(node, ts, Sm_t1, Sf1_t1, Sf2_t1, Sf3_t1, Ss1_t1, Q_t)

    return Q_t
end

"""
    update_params!(node::HyModNode, Sm_max::F, B::F, alpha::F, Kf::F, Ks::F) where {F<:Float64}

Update parameters for HyMod.
"""
function update_params!(node::HyModNode, Sm_max::F, B::F, alpha::F, Kf::F, Ks::F) where {F<:Float64}
    node.Sm_max = Param(Sm_max, bounds=node.Sm_max.bounds::Tuple)
    node.B = Param(B, bounds=node.B.bounds::Tuple)
    node.alpha = Param(alpha, bounds=node.alpha.bounds::Tuple)
    node.Kf = Param(Kf, bounds=node.Kf.bounds::Tuple)
    node.Ks = Param(Ks, bounds=node.Ks.bounds::Tuple)
end

function update_state!(
    node::HyModNode, ts::Int, Sm::F, Sf1::F, Sf2::F, Sf3::F, Ss1::F, Q::F
)::Nothing where {F<:Float64}
    # State for time t+1
    node.Sm[ts+1] = Sm
    node.Sf1[ts+1] = Sf1
    node.Sf2[ts+1] = Sf2
    node.Sf3[ts+1] = Sf3
    node.Ss1[ts+1] = Ss1

    # Outflow for time t
    node.outflow[ts] = Q

    return nothing
end

function reset!(node::HyModNode)::Nothing
    node.Sm = Float64[node.Sm[1]]
    node.Sf1 = Float64[node.Sf1[1]]
    node.Sf2 = Float64[node.Sf2[1]]
    node.Sf3 = Float64[node.Sf3[1]]
    node.Ss1 = Float64[node.Ss1[1]]
    node.outflow = Float64[]

    return nothing
end
