using Parameters
using ModelParameters

include("./components/p_and_et.jl")
include("./components/cmd.jl")
include("./components/flow.jl")

abstract type IHACRESNode <: NetworkNode end


Base.@kwdef mutable struct IHACRESBilinearNode{P, A<:AbstractFloat} <: IHACRESNode
    const name::String
    const area::A

    # https://wiki.ewater.org.au/display/SD41/IHACRES-CMD+-+SRG
    d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    e::P = Param(1.0, bounds=(0.999, 1.0))  # temperature to PET conversion factor
    f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    a::P = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == (1/tau_q)
    b::P = Param(0.1, bounds=(1e-3, 0.1))  # slowflow storage coefficent == (1/tau_s)

    storage_coef::P = Param(2.9, bounds=(1e-10, 10.0))
    alpha::P = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    # const level_params::Array{P, 1} = [
    #     Param(-0.01, bounds=(-10.0, -0.01)),  # p1
    #     Param(0.8, bounds=(0.0, 1.5)),  # p2
    #     Param(4.5, bounds=(0.0, 20.0)), # p3
    #     Param(5.0, bounds=(1.0, 10.0)), # p4
    #     Param(0.35, bounds=(0.0, 1.0)), # p5
    #     Param(1.41, bounds=(-2.0, 2.0)), # p6
    #     Param(-1.45, bounds=(-2.5, 0.0)), # p7
    #     Param(6.75, bounds=(0.0, 10.0)), # p8
    #     Param(150.0, bounds=(50.0, 200.0)) # CTF (height in local datum where stream will cease to flow)
    # ]

    storage::Array{A} = [100.0]  # CMD
    quick_store::Array{A} = [0.0]
    slow_store::Array{A} = [0.0]
    outflow::Array{A} = []
    effective_rainfall::Array{A} = []
    et::Array{A} = []
    inflow::Array{A} = []
    # level::Array{A} = []
    gw_store::Array{A} = [0.0]

    obj_func::Function = obj_func
end


function prep_state!(node::IHACRESNode, timesteps::Int64)::Nothing
    resize!(node.storage, timesteps+1)
    node.storage[2:end] .= 0.0

    resize!(node.quick_store, timesteps+1)
    node.quick_store[2:end] .= 0.0

    resize!(node.slow_store, timesteps+1)
    node.slow_store[2:end] .= 0.0

    node.outflow = fill(0.0, timesteps)
    node.effective_rainfall = fill(0.0, timesteps)
    node.et = fill(0.0, timesteps)
    node.inflow = fill(0.0, timesteps)
    # node.level = fill(0.0, timesteps)

    resize!(node.gw_store, timesteps+1)
    node.gw_store[2:end] .= 0.0

    return nothing
end


function IHACRESBilinearNode(name::String, spec::AbstractDict)
    n = create_node(IHACRESBilinearNode, name, spec["area"])
    node_params = copy(spec["parameters"])

    # Initial Catchment Moisture Deficit
    n.storage = [spec["initial_storage"]]

    # if haskey(spec, "level_params")
    #     n_lparams = n.level_params
    #     s_lparams = spec["level_params"]
    #     node_params["level_params"] = Param[
    #         Param(s_lparams[1], bounds=n_lparams[1].bounds)
    #         Param(s_lparams[2], bounds=n_lparams[2].bounds)
    #         Param(s_lparams[3], bounds=n_lparams[3].bounds)
    #         Param(s_lparams[4], bounds=n_lparams[4].bounds)
    #         Param(s_lparams[5], bounds=n_lparams[5].bounds)
    #         Param(s_lparams[6], bounds=n_lparams[6].bounds)
    #         Param(s_lparams[7], bounds=n_lparams[7].bounds)
    #         Param(s_lparams[8], bounds=n_lparams[8].bounds)
    #         Param(s_lparams[9], bounds=n_lparams[9].bounds)
    #     ]
    # end

    for (k, p) in node_params
        s = Symbol(k)
        if p isa String
            p = eval(Meta.parse(p))
        end

        try
            f = getfield(n, s)
            setfield!(n, s, Param(p, bounds=f.bounds))
        catch err
            msg = sprint(showerror, err, catch_backtrace())
            if occursin("no field bounds", string(msg))
                setfield!(n, s, p)
            else
                throw(err)
            end
        end
    end

    return n
end


"""
    IHACRESBilinearNode(name::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                store::Float64, quick::Float64, slow::Float64, gw_store::Float64)

Create a IHACRES node that adopts the bilinear CMD module.
"""
function IHACRESBilinearNode(name::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                      a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                      store::Float64, quick::Float64, slow::Float64, gw_store::Float64)
    n = create_node(IHACRESBilinearNode, name, area)

    update_params!(n, d, d2, e, f, a, b, s_coef, alpha)

    n.storage = [store]
    n.quick_store = [quick]
    n.slow_store = [slow]
    n.gw_store = [gw_store]

    return n
end


function update_state!(
    s_node::IHACRESNode,
    storage::F,
    e_rainfall::F,
    et::F,
    qflow_store::F,
    sflow_store::F,
    outflow::F,
    gw_store::F
)::Nothing where {F<:Float64}
    push!(s_node.storage, storage)
    push!(s_node.effective_rainfall, e_rainfall)
    push!(s_node.et, et)
    push!(s_node.outflow, outflow)

    push!(s_node.quick_store, qflow_store)
    push!(s_node.slow_store, sflow_store)
    # push!(s_node.level, level)
    push!(s_node.gw_store, gw_store)

    return nothing
end
function update_state!(
    s_node::IHACRESNode,
    ts::Int64,
    storage::Float64,
    e_rainfall::Float64,
    et::Float64,
    qflow_store::Float64,
    sflow_store::Float64,
    inflow::Float64,
    outflow::Float64,
    gw_store::Float64
)::Nothing
    s_node.storage[ts+1] = storage
    s_node.quick_store[ts+1] = qflow_store
    s_node.slow_store[ts+1] = sflow_store
    s_node.gw_store[ts+1] = gw_store

    s_node.inflow[ts] = inflow
    s_node.outflow[ts] = outflow

    s_node.effective_rainfall[ts] = e_rainfall
    s_node.et[ts] = et
    # s_node.level[ts] = level

    return nothing
end

"""
    run_timestep!(
        node::IHACRESBilinearNode, climate::Climate, timestep::Int,
        inflow::Float64, extraction::Float64, exchange::Float64
    )

Run the given IHACRESBilinearNode for a timestep.
"""
function run_timestep!(
    node::IHACRESNode, climate::Climate, timestep::Int64;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Float64
    ts = timestep
    P, ET = climate_values(node, climate, ts)

    return run_timestep!(node, ts, P, ET, in_flow, ext, flux)
end
function run_timestep!(
    node::IHACRESNode, P::F, ET::F, ts::Int64;
    inflow=nothing, extraction=nothing, exchange=nothing
)::F where {F<:AbstractFloat}
    in_flow = timestep_value(ts, node.name, "inflow", inflow)
    ext = timestep_value(ts, node.name, "extraction", extraction)
    flux = timestep_value(ts, node.name, "exchange", exchange)

    return run_ihacres!(node, P, ET, ts, in_flow, ext, flux)
end

function run_ihacres!(
    s_node::IHACRESNode,
    rain::F,
    evap::F,
    ts::Int64,
    inflow::F,
    ext::F,
    gw_exchange::F
)::F where {F<:Float64}
    current_store::F = s_node.storage[ts]
    quick_store::F = s_node.quick_store[ts]
    slow_store::F = s_node.slow_store[ts]
    gw_store::F = s_node.gw_store[ts]

    (mf, e_rainfall, recharge) = calc_ft_interim_cmd(
        current_store,
        rain,
        s_node.d.val::F,
        s_node.d2.val::F,
        s_node.alpha.val::F
    )

    et::F = calc_ET_from_E(
        s_node.e.val::F,
        evap,
        mf,
        s_node.f.val::F,
        s_node.d.val::F
    )

    cmd::F = calc_cmd(
        current_store,
        rain,
        et,
        e_rainfall,
        recharge
    )

    (nq_store, ns_store, outflow) = calc_ft_flows(
        quick_store,
        slow_store,
        e_rainfall,
        recharge,
        s_node.area,
        s_node.a.val::F,
        s_node.b.val::F
    )

    (gw_store, outflow) = routing(
        gw_store,
        s_node.storage_coef.val::F,
        inflow,
        outflow,
        ext,
        gw_exchange
    )

    # level::Float64 = calc_ft_level(
    #     outflow,
    #     Vector{Float64}(s_node.level_params)
    # )

    update_state!(s_node, ts, cmd, e_rainfall, et, nq_store, ns_store, inflow, outflow, gw_store)

    return outflow
end


"""
    param_info(node::IHACRESNode)::Tuple

Extract node parameter names, values, and bounds for IHACRESNode types.
"""
function param_info(node::IHACRESNode)::Tuple
    tmp = Model(node)
    values = collect(tmp[:val])
    bounds = collect(tmp[:bounds])
    param_names = collect(tmp[:fieldname])

    # if with_level
    #     level_param_vals = map(x->x.val, node.level_params)
    #     append!(values, level_param_vals)

    #     level_param_bounds = map(x->x.bounds, node.level_params)
    #     append!(bounds, level_param_bounds)
    # end

    return param_names, values, bounds
end


function run_node_with_temp!(sn::StreamfallNetwork, nid::Int64, climate::Climate;
                             inflow=nothing, extraction=nothing, exchange=nothing)
    node = sn[nid]
    timesteps = sim_length(climate)
    @inbounds for ts in 1:timesteps
        run_node_with_temp!(node, climate, ts;
                            inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow, node.level
end


function run_node_with_temp!(node::IHACRESBilinearNode, climate::Climate;
                             inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    @inbounds for ts in 1:timesteps
        run_node_with_temp!(node, climate, ts;
                            inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow, node.level
end


function run_node_with_temp!(node::IHACRESBilinearNode, climate::Climate, timestep::Int64;
                             inflow=nothing, extraction=nothing, exchange=nothing)
    ts = timestep
    P, T = climate_values(node, climate, ts)
    i = timestep_value(ts, node.name, "_inflow", inflow)
    ext = timestep_value(ts, node.name, "_extraction", extraction)
    flux = timestep_value(ts, node.name, "_exchange", exchange)
    run_node_with_temp!(node, P, T, i, ext, flux)
end


"""
    run_node_with_temp!(s_node::IHACRESBilinearNode,
                        rain::Float64,
                        temp::Float64,
                        inflow::Float64,
                        ext::Float64,
                        gw_exchange::Float64=0.0;
                        current_store=nothing,
                        quick_store=nothing,
                        slow_store=nothing
                        gw_store=nothing)::Tuple{Float64, Float64}

Run node with temperature data to calculate outflow and update state.
"""
function run_node_with_temp!(s_node::IHACRESBilinearNode,
                             rain::Float64,
                             temp::Float64,
                             inflow::Float64,
                             ext::Float64,
                             gw_exchange::Float64=0.0;
                             current_store=nothing,
                             quick_store=nothing,
                             slow_store=nothing,
                             gw_store=nothing)::Tuple{Float64, Float64}
    current_store = s_node.storage[end]
    quick_store = s_node.quick_store[end]
    slow_store = s_node.slow_store[end]
    gw_store = s_node.gw_store[end]

    (mf, e_rainfall, recharge) = calc_ft_interim_cmd(
        interim_results,
        current_store,
        rain,
        s_node.d,
        s_node.d2,
        s_node.alpha
    )

    et::Float64 = calc_ET_from_T(
        s_node.e,
        temp,
        mf,
        s_node.f,
        s_node.d
    )

    cmd::Float64 = calc_cmd(
        mf,
        rain,
        et,
        e_rainfall,
        recharge
    )

    (nq_store, ns_store, outflow) = calc_ft_flows(
        quick_store,
        slow_store,
        e_rainfall,
        recharge,
        s_node.area,
        s_node.a,
        s_node.b
    )

    @assert any(isnan.((nq_store, ns_store, outflow))) == false

    (gw_store, outflow) = routing(
        gw_store,
        s_node.storage_coef,
        inflow,
        outflow,
        ext,
        gw_exchange
    )

    # level_params = Array{Float64}(s_node.level_params)
    # level::Float64 = calc_ft_level(
    #     outflow,
    #     level_params
    # )

    push!(s_node.inflow, inflow)
    update_state(s_node, cmd, e_rainfall, et, nq_store, ns_store, outflow, gw_store)

    return outflow, level
end


"""
    update_params!(
        node::IHACRESBilinearNode,
        d::F,
        d2::F,
        e::F,
        f::F,
        a::F,
        b::F,
        s_coef::F,
        alpha::F
    )::Nothing where {F<:Float64}

Update model parameters.
"""
function update_params!(
    node::IHACRESBilinearNode,
    d::F,
    d2::F,
    e::F,
    f::F,
    a::F,
    b::F,
    s_coef::F,
    alpha::F
)::Nothing where {F<:Float64}
    node.d = Param(d, bounds=node.d.bounds::Tuple)
    node.d2 = Param(d2, bounds=node.d2.bounds::Tuple)
    node.e = Param(e, bounds=node.e.bounds::Tuple)
    node.f = Param(f, bounds=node.f.bounds::Tuple)
    node.a = Param(a, bounds=node.a.bounds::Tuple)
    node.b = Param(b, bounds=node.b.bounds::Tuple)
    node.storage_coef = Param(s_coef, bounds=node.storage_coef.bounds::Tuple)
    node.alpha = Param(alpha, bounds=node.alpha.bounds::Tuple)

    return nothing
end

"""
    reset!(s_node::IHACRESNode)::Nothing

Reset node. Clears all states back to their initial values.
"""
function reset!(s_node::IHACRESNode)::Nothing
    s_node.storage = [s_node.storage[1]]
    s_node.quick_store = [s_node.quick_store[1]]
    s_node.slow_store = [s_node.slow_store[1]]
    s_node.gw_store = [s_node.gw_store[1]]
    s_node.outflow = []
    # s_node.level = []
    s_node.effective_rainfall = []
    s_node.et = []
    s_node.inflow = []

    return nothing
end

"""
    extract_spec!(node::DamNode, spec::AbstractDict)::Nothing

Additional processing to extract IHACRES-specific details.
"""
function extract_spec!(node::IHACRESNode, spec::AbstractDict)::Nothing
    spec["initial_storage"] = node.storage[1]

    return nothing
end
