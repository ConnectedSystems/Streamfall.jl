using Parameters
using ModelParameters


Base.@kwdef mutable struct ExpuhNode{P, A<:AbstractFloat} <: IHACRESNode
    name::String
    area::A

    # https://wiki.ewater.org.au/display/SD41/IHACRES-CMD+-+SRG
    d::P = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    d2::P = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold, multiplier applied to d
    e::P = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    f::P = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    tau_q::P = Param(1.0, bounds=(0.0, 5.0))
    tau_s::P = Param(5.0, bounds=(5.0, 200.0))
    v_s::P = Param(0.5, bounds=(0.0, 1.0))

    storage_coef::P = Param(2.9, bounds=(0.2, 10.0))

    # level_params::Array{P, 1} = [
    #     Param(-0.01, bounds=(-10.0, -0.01)),  # p1
    #     Param(0.8, bounds=(0.0, 1.5)),  # p2
    #     Param(4.5, bounds=(0.0, 20.0)), # p3
    #     Param(5.0, bounds=(1.0, 10.0)), # p4
    #     Param(0.35, bounds=(0.0, 1.0)), # p5
    #     Param(1.41, bounds=(-2.0, 2.0)), # p6
    #     Param(-1.45, bounds=(-2.5, 0.0)), # p7
    #     Param(6.75, bounds=(0.0, 10.0)), # p8
    #     Param(150.0, bounds=(50.0, 200.0)) # ctf
    # ]

    storage::Array{A} = [100.0]
    quick_store::Array{A} = [0.0]
    slow_store::Array{A} = [0.0]
    outflow::Array{A} = []
    effective_rainfall::Array{A} = []
    et::Array{A} = []
    inflow::Array{A} = []
    level::Array{A} = []
    gw_store::Array{A} = [0.0]

    obj_func::Function = obj_func
end


function ExpuhNode(name::String, spec::Dict)
    n = create_node(ExpuhNode, name, spec["area"])

    node_params = spec["parameters"]
    n_lparams = n.level_params
    s_lparams = spec["level_params"]
    # node_params["level_params"] = Param[
    #     Param(s_lparams[1], bounds=n_lparams[1].bounds)
    #     Param(s_lparams[2], bounds=n_lparams[2].bounds)
    #     Param(s_lparams[3], bounds=n_lparams[3].bounds)
    #     Param(s_lparams[4], bounds=n_lparams[4].bounds)
    #     Param(s_lparams[5], bounds=n_lparams[5].bounds)
    #     Param(s_lparams[6], bounds=n_lparams[6].bounds)
    #     Param(s_lparams[7], bounds=n_lparams[7].bounds)
    #     Param(s_lparams[8], bounds=n_lparams[8].bounds)
    #     Param(s_lparams[9], bounds=n_lparams[9].bounds)
    # ]

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
            if occursin("no field bounds", msg)
                setfield!(n, s, p)
            else
                throw(err)
            end
        end
    end

    return n
end


function ExpuhNode(name::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                    tau_q::Float64, tau_s::Float64, v_s::Float64, s_coef::Float64,
                    store::Float64, quick::Float64, slow::Float64)
    return ExpuhNode{Param, Float64}(
        name=name,
        area=area,
        d=d,
        d2=d2,
        e=e,
        f=f,
        tau_q=tau_q,
        tau_s=tau_s,
        v_s=v_s,
        storage_coef=s_coef,
        storage=[store],
        quick_store=[quick],
        slow_store=[slow]
    )
end

"""
    run_node!(s_node::ExpuhNode, rain::Float64, evap::Float64,
              inflow::Float64, ext::Float64, gw_exchange::Float64;
              current_store=nothing,
              quick_store=nothing,
              slow_store=nothing,
            )::Tuple

Run given IHACRES ExpuhNode for a time step based on last known state.
"""
function run_node!(s_node::ExpuhNode,
                   rain::Float64,
                   evap::Float64,
                   inflow::Float64,
                   ext::Float64,
                   gw_exchange::Float64;
                   current_store=nothing,
                   quick_store=nothing,
                   slow_store=nothing,
                )::Tuple

    if !isnothing(current_store)
        s_node.storage[end] = current_store
    end

    if !isnothing(quick_store)
        s_node.quick_store[end] = quick_store
    end

    if !isnothing(slow_store)
        s_node.slow_store[end] = slow_store
    end

    cmd::Float64 = s_node.storage[end]

    e_rainfall = @ccall IHACRES.calc_effective_rainfall(rain::Cdouble, cmd::Cdouble, s_node.d::Cdouble, s_node.d2::Cdouble)::Cdouble
    mf = @ccall IHACRES.calc_trig_interim_cmd(cmd::Cdouble, s_node.d::Cdouble, e_rainfall::Cdouble)::Cdouble

    et::Float64 = @ccall IHACRES.calc_ET(s_node.e::Cdouble, evap::Cdouble, mf::Cdouble, s_node.f::Cdouble, s_node.d::Cdouble)::Cdouble
    cmd = @ccall IHACRES.calc_cmd(mf::Cdouble, rain::Cdouble, et::Cdouble, e_rainfall::Cdouble)::Cdouble

    (prev_q, prev_s) = (s_node.quick_store[end], s_node.slow_store[end])

    flow_res = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_flows(flow_res::Ptr{Cdouble}, prev_q::Cdouble, prev_s::Cdouble, s_node.v_s::Cdouble, e_rainfall::Cdouble,
                              s_node.area::Cdouble, s_node.tau_q::Cdouble, s_node.tau_s::Cdouble)::Cvoid
    (quick_store, slow_store, outflow) = flow_res

    gw_store = s_node.gw_store[end]
    routing_res = [0.0, 0.0]
    @ccall IHACRES.routing(
            routing_res::Ptr{Cdouble},
            gw_store::Cdouble,
            s_node.storage_coef::Cdouble,
            inflow::Cdouble,
            outflow::Cdouble,
            ext::Cdouble,
            gw_exchange::Cdouble)::Cvoid
    (gw_store, outflow) = routing_res

    level::Float64 = @ccall IHACRES.calc_ft_level(outflow::Cdouble, s_node.level_params::Ptr{Cdouble})::Cdouble

    update_state!(s_node, cmd, e_rainfall, et, quick_store, slow_store, outflow, level, gw_store)

    return (outflow, level)
end


function run_node_with_temp!(s_node::ExpuhNode, rain::Float64, temp::Float64, inflow::Float64, ext::Float64; gw_exchange::Float64=0.0)::Tuple

    cmd::Float64 = s_node.storage[end]

    # rain = rain * s_node.area
    e_rainfall = @ccall IHACRES.calc_effective_rainfall(rain::Cdouble, cmd::Cdouble, s_node.d::Cdouble, s_node.d2::Cdouble)::Cdouble
    # mf = @ccall IHACRES.calc_linear_interim_cmd(cmd::Cdouble, s_node.d::Cdouble, e_rainfall::Cdouble)::Cdouble
    mf = @ccall IHACRES.calc_trig_interim_cmd(cmd::Cdouble, s_node.d::Cdouble, e_rainfall::Cdouble)::Cdouble

    et::Float64 = @ccall IHACRES.calc_ET_from_T(
        s_node.e::Cdouble,
        temp::Cdouble,
        mf::Cdouble,
        s_node.f::Cdouble,
        s_node.d::Cdouble
    )::Cdouble

    cmd = @ccall IHACRES.calc_cmd(current_store::Cdouble, rain::Cdouble, et::Cdouble, e_rainfall::Cdouble)::Cdouble

    (prev_q, prev_s) = (s_node.quick_store[end], s_node.slow_store[end])

    flow_res = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_flows(flow_res::Ptr{Cdouble}, prev_q::Cdouble, prev_s::Cdouble, s_node.v_s::Cdouble, e_rainfall::Cdouble,
                              s_node.area::Cdouble, s_node.tau_q::Cdouble, s_node.tau_s::Cdouble)::Cvoid
    (quick_store, slow_store, outflow) = flow_res

    gw_store = s_node.gw_store[end]
    routing_res = [0.0, 0.0]
    @ccall IHACRES.routing(
            routing_res::Ptr{Cdouble},
            gw_store::Cdouble,
            s_node.storage_coef::Cdouble,
            inflow::Cdouble,
            outflow::Cdouble,
            ext::Cdouble,
            gw_exchange::Cdouble)::Cvoid
    (gw_store, outflow) = routing_res

    level::Float64 = @ccall IHACRES.calc_ft_level(outflow::Cdouble, s_node.level_params::Ptr{Cdouble})::Cdouble

    update_state!(s_node, cmd, e_rainfall, et, quick_store, slow_store, outflow, level, gw_store)

    return (outflow, level)
end


function update_params!(node::ExpuhNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                        tau_q::Float64, tau_s::Float64, v_s::Float64, s_coef::Float64)::Nothing
    node.d = Param(d, bounds=node.d.bounds)
    node.d2 = Param(d2, bounds=node.d2.bounds)
    node.e = Param(e, bounds=node.e.bounds)
    node.f = Param(f, bounds=node.f.bounds)
    node.tau_q = Param(tau_q, bounds=node.tau_q.bounds)
    node.tau_s = Param(tau_s, bounds=node.tau_s.bounds)
    node.v_s = Param(v_s, bounds=node.v_s.bounds)
    node.storage_coef = Param(s_coef, bounds=node.storage_coef.bounds)

    return nothing
end
