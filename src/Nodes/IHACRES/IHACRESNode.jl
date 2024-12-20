using Parameters
using ModelParameters


abstract type IHACRESNode <: NetworkNode end


Base.@kwdef mutable struct BilinearNode{P, A<:AbstractFloat} <: IHACRESNode
    const name::String
    const area::A

    # https://wiki.ewater.org.au/display/SD41/IHACRES-CMD+-+SRG
    d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    e::P = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    a::P = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == (1/tau_q)
    b::P = Param(0.1, bounds=(1e-3, 0.1))  # slowflow storage coefficent == (1/tau_s)

    # These changes require ihacres_nim to be changed as well...
    # a::A = Param(0.9, bounds=(0.1, 10.0))  # quickflow storage coefficient == exp(-1/tau_q)
    # b::A = Param(0.1, bounds=(10.0, 1000.0))  # slowflow storage coefficent == exp(-1/tau_s)
    # storage_threshold::A = Param(10.0, ...)  # optional threshold controlling bounds between `a` and `b`
    storage_coef::P = Param(2.9, bounds=(1e-10, 10.0))
    alpha::P = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    const level_params::Array{P, 1} = [
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
    node.level = fill(0.0, timesteps)

    resize!(node.gw_store, timesteps+1)
    node.gw_store[2:end] .= 0.0

    return nothing
end


function BilinearNode(name::String, spec::Dict)
    n = create_node(BilinearNode, name, spec["area"])
    node_params = copy(spec["parameters"])

    if haskey(spec, "level_params")
        n_lparams = n.level_params
        s_lparams = spec["level_params"]
        node_params["level_params"] = Param[
            Param(s_lparams[1], bounds=n_lparams[1].bounds)
            Param(s_lparams[2], bounds=n_lparams[2].bounds)
            Param(s_lparams[3], bounds=n_lparams[3].bounds)
            Param(s_lparams[4], bounds=n_lparams[4].bounds)
            Param(s_lparams[5], bounds=n_lparams[5].bounds)
            Param(s_lparams[6], bounds=n_lparams[6].bounds)
            Param(s_lparams[7], bounds=n_lparams[7].bounds)
            Param(s_lparams[8], bounds=n_lparams[8].bounds)
            Param(s_lparams[9], bounds=n_lparams[9].bounds)
        ]
    end

    node_params["storage"] = [node_params["initial_storage"]]
    delete!(node_params, "initial_storage")

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
    BilinearNode(name::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                store::Float64, quick::Float64, slow::Float64, gw_store::Float64)

Create a IHACRES node that adopts the bilinear CMD module.
"""
function BilinearNode(name::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                      a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                      store::Float64, quick::Float64, slow::Float64, gw_store::Float64)
    n = create_node(BilinearNode, name, area)

    update_params!(n, d, d2, e, f, a, b, s_coef, alpha)

    n.storage = [store]
    n.quick_store = [quick]
    n.slow_store = [slow]
    n.gw_store = [gw_store]

    return n
end


function update_state!(s_node::IHACRESNode, storage::Float64, e_rainfall::Float64, et::Float64,
                      qflow_store::Float64, sflow_store::Float64, outflow::Float64,
                      level::Float64, gw_store::Float64)::Nothing
    push!(s_node.storage, storage)
    push!(s_node.effective_rainfall, e_rainfall)
    push!(s_node.et, et)
    push!(s_node.outflow, outflow)

    push!(s_node.quick_store, qflow_store)
    push!(s_node.slow_store, sflow_store)
    push!(s_node.level, level)
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
    level::Float64,
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
    s_node.level[ts] = level

    return nothing
end

"""
    run_node!(node::IHACRESNode, climate::Climate, timestep::Int;
                  inflow=nothing, extraction=nothing, exchange=nothing)

Run a specific node for a specified time step.

# Arguments
- `node::IHACRESNode` :
- `climate::Climate` :
- `ts::Int` : current time step
- `inflow::DataFrame` : Time series of inflows from any upstream node.
- `extraction::DataFrame` : Time series of water orders (expects column of `_releases`)
- `exchange::DataFrame` : Time series of groundwater flux
"""
function run_node!(node::IHACRESNode, climate::Climate, ts::Int;
                   inflow=nothing, extraction=nothing, exchange=nothing)
    if checkbounds(Bool, node.outflow, ts)
        if node.outflow[ts] != undef
            # already ran for this time step so no need to run
            return node.outflow[ts], node.level[ts]
        end
    end

    node_name = node.name
    rain, et = climate_values(node, climate, ts)
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)

    return run_timestep!(node, rain, et, ts; inflow=in_flow, extraction=wo, exchange=ex)
end

"""
    run_timestep!(
        node::BilinearNode, climate::Climate, timestep::Int,
        inflow::Float64, extraction::Float64, exchange::Float64
    )

Run the given BilinearNode for a time step.
"""
function run_timestep!(
    node::IHACRESNode, climate::Climate, timestep::Int;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Float64
    ts = timestep
    P, ET = climate_values(node, climate, ts)

    return run_timestep!(node, ts, P, ET, in_flow, ext, flux)
end
function run_timestep!(node, P, ET, ts; inflow=nothing, extraction=nothing, exchange=nothing)
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
)::Tuple{F, F} where {F<:Float64}
    current_store::F = s_node.storage[ts]
    quick_store::F = s_node.quick_store[ts]
    slow_store::F = s_node.slow_store[ts]
    gw_store::F = s_node.gw_store[ts]

    cache_res = s_node.cache.i_cache
    @ccall IHACRES.calc_ft_interim_cmd(
        cache_res::Ptr{Cdouble},
        current_store::Cdouble,
        rain::Cdouble,
        s_node.d.val::Cdouble,
        s_node.d2.val::Cdouble,
        s_node.alpha.val::Cdouble
    )::Cvoid
    (mf, e_rainfall, recharge) = cache_res

    et::F = @ccall IHACRES.calc_ET(
        s_node.e.val::Cdouble,
        evap::Cdouble,
        mf::Cdouble,
        s_node.f.val::Cdouble,
        s_node.d.val::Cdouble
    )::Cdouble

    cmd::F = @ccall IHACRES.calc_cmd(
        current_store::Cdouble,
        rain::Cdouble,
        et::Cdouble,
        e_rainfall::Cdouble,
        recharge::Cdouble
    )::Cdouble

    @ccall IHACRES.calc_ft_flows(
        cache_res::Ptr{Cdouble},
        quick_store::Cdouble,
        slow_store::Cdouble,
        e_rainfall::Cdouble,
        recharge::Cdouble,
        s_node.area::Cdouble,
        s_node.a.val::Cdouble,
        s_node.b.val::Cdouble
    )::Cvoid
    (nq_store, ns_store, outflow) = cache_res

    routing_res = s_node.cache.r_cache
    @ccall IHACRES.routing(
        routing_res::Ptr{Cdouble},
        gw_store::Cdouble,
        s_node.storage_coef.val::Cdouble,
        inflow::Cdouble,
        outflow::Cdouble,
        ext::Cdouble,
        gw_exchange::Cdouble
    )::Cvoid

    (gw_store, outflow) = routing_res

    level_params = Vector{Float64}(s_node.level_params)
    level::Float64 = @ccall IHACRES.calc_ft_level(
        outflow::Cdouble,
        level_params::Ptr{Cdouble}
    )::Cdouble

    update_state!(s_node, ts, cmd, e_rainfall, et, nq_store, ns_store, inflow, outflow, level, gw_store)

    return outflow, level
end


"""
    param_info(node::IHACRESNode; with_level::Bool = true)::Tuple

Extract node parameter names, values, and bounds for IHACRESNode types.
"""
function param_info(node::IHACRESNode; with_level::Bool = false)::Tuple
    tmp = Model(node)
    values = collect(tmp[:val])
    bounds = collect(tmp[:bounds])
    param_names = collect(tmp[:fieldname])

    if with_level
        level_param_vals = map(x->x.val, node.level_params)
        append!(values, level_param_vals)

        level_param_bounds = map(x->x.bounds, node.level_params)
        append!(bounds, level_param_bounds)
    end

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


function run_node_with_temp!(node::BilinearNode, climate::Climate;
                             inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    @inbounds for ts in 1:timesteps
        run_node_with_temp!(node, climate, ts;
                            inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow, node.level
end


function run_node_with_temp!(node::BilinearNode, climate::Climate, timestep::Int64;
                             inflow=nothing, extraction=nothing, exchange=nothing)
    ts = timestep
    P, T = climate_values(node, climate, ts)
    i = timestep_value(ts, node.name, "_inflow", inflow)
    ext = timestep_value(ts, node.name, "_extraction", extraction)
    flux = timestep_value(ts, node.name, "_exchange", exchange)
    run_node_with_temp!(node, P, T, i, ext, flux)
end


"""
    run_node_with_temp!(s_node::BilinearNode,
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
function run_node_with_temp!(s_node::BilinearNode,
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

    interim_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_interim_cmd(interim_results::Ptr{Cdouble},
                                       current_store::Cdouble,
                                       rain::Cdouble,
                                       s_node.d::Cdouble,
                                       s_node.d2::Cdouble,
                                       s_node.alpha::Cdouble)::Cvoid

    @assert any(isnan.(interim_results)) == false
    (mf, e_rainfall, recharge) = interim_results

    et::Float64 = @ccall IHACRES.calc_ET_from_T(
        s_node.e::Cdouble,
        temp::Cdouble,
        mf::Cdouble,
        s_node.f::Cdouble,
        s_node.d::Cdouble
    )::Cdouble

    cmd::Float64 = @ccall IHACRES.calc_cmd(
        mf::Cdouble,
        rain::Cdouble,
        et::Cdouble,
        e_rainfall::Cdouble,
        recharge::Cdouble
    )::Cdouble

    flow_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_flows(
        flow_results::Ptr{Cdouble},
        quick_store::Cdouble,
        slow_store::Cdouble,
        e_rainfall::Cdouble,
        recharge::Cdouble,
        s_node.area::Cdouble,
        s_node.a::Cdouble,
        s_node.b::Cdouble
    )::Cvoid

    @assert any(isnan.(flow_results)) == false
    (nq_store, ns_store, outflow) = flow_results

    routing_res = [0.0, 0.0]
    @ccall IHACRES.routing(
        routing_res::Ptr{Cdouble},
        gw_store::Cdouble,
        s_node.storage_coef::Cdouble,
        inflow::Cdouble,
        outflow::Cdouble,
        ext::Cdouble,
        gw_exchange::Cdouble)::Cvoid

    @assert any(isnan.(routing_res)) == false
    (gw_store, outflow) = routing_res

    level_params = Array{Float64}(s_node.level_params)
    level::Float64 = @ccall IHACRES.calc_ft_level(
        outflow::Cdouble,
        level_params::Ptr{Cdouble}
    )::Cdouble

    push!(s_node.inflow, inflow)
    update_state(s_node, cmd, e_rainfall, et, nq_store, ns_store, outflow, level, gw_store)

    return outflow, level
end


"""
    update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                   a::Float64, b::Float64, s_coef::Float64, alpha::Float64)::Nothing

Update model parameters.
"""
function update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                        a::Float64, b::Float64, s_coef::Float64, alpha::Float64)::Nothing
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
    update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                   a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                   p1::Float64, p2::Float64, p3::Float64, p4::Float64, p5::Float64, p6::Float64, p7::Float64, p8::Float64, CTF::Float64)::Nothing

Update all parameters.
"""
function update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                        a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                        p1::Float64, p2::Float64, p3::Float64, p4::Float64, p5::Float64,
                        p6::Float64, p7::Float64, p8::Float64, CTF::Float64)::Nothing
    node.d = Param(d, bounds=node.d.bounds::Tuple)
    node.d2 = Param(d2, bounds=node.d2.bounds::Tuple)
    node.e = Param(e, bounds=node.e.bounds::Tuple)
    node.f = Param(f, bounds=node.f.bounds::Tuple)
    node.a = Param(a, bounds=node.a.bounds::Tuple)
    node.b = Param(b, bounds=node.b.bounds::Tuple)
    node.storage_coef = Param(s_coef, bounds=node.storage_coef.bounds::Tuple)
    node.alpha = Param(alpha, bounds=node.alpha.bounds::Tuple)

    n_lparams = node.level_params
    node.level_params = [
        Param(p1, bounds=n_lparams[1].bounds::Tuple)
        Param(p2, bounds=n_lparams[2].bounds::Tuple)
        Param(p3, bounds=n_lparams[3].bounds::Tuple)
        Param(p4, bounds=n_lparams[4].bounds::Tuple)
        Param(p5, bounds=n_lparams[5].bounds::Tuple)
        Param(p6, bounds=n_lparams[6].bounds::Tuple)
        Param(p7, bounds=n_lparams[7].bounds::Tuple)
        Param(p8, bounds=n_lparams[8].bounds::Tuple)
        Param(CTF, bounds=n_lparams[9].bounds::Tuple)
    ]

    return nothing
end

# TODO
# update_bounds!()


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
    s_node.level = []
    s_node.effective_rainfall = []
    s_node.et = []
    s_node.inflow = []

    return nothing
end