using Parameters
using ModelParameters


Base.@kwdef mutable struct IHACRESNode{A} <: NetworkNode{A}
    @network_node

    # https://wiki.ewater.org.au/display/SD41/IHACRES-CMD+-+SRG
    d::A = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    d2::A = Param(2.0, bounds=(10.0, 500.0))   # flow threshold2
    e::A = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    f::A = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    a::A = Param(200.0, bounds=(0.1, 100.0))
    b::A = Param(200.0, bounds=(0.001, 1.0))
    storage_coef::A = Param(2.9, bounds=(0.2, 10.0))
    alpha::A = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    level_params::Array{A} = [
        Param(-0.01, bounds=(-0.01, -10.0))
        Param(0.8, bounds=(0.0, 1.5))
        Param(4.5, bounds=(0.0, 20.0))
        Param(5.0, bounds=(1.0, 10.0))
        Param(0.35, bounds=(0.0, 1.0))
        Param(1.41, bounds=(-2.0, 2.0))
        Param(-1.45, bounds=(0.0, -2.5))
        Param(6.75, bounds=(0.0, 10.0))
        Param(150.0, bounds=(50.0, 200.0))
    ]

    storage::Array{Float64} = [100.0]
    quickflow::Array{Float64} = [100.0]
    slowflow::Array{Float64} = [100.0]
    outflow::Array{Float64} = []
    effective_rainfall::Array{Float64} = []
    et::Array{Float64} = []
    inflow::Array{Float64} = []
    level::Array{Float64} = []
end


function IHACRESNode(node_id::String, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64, 
                    a::Float64, b::Float64, s_coef::Float64, alpha::Float64, 
                    store::Float64, quick::Float64, slow::Float64)
    return IHACRESNode(
        node_id=node_id,
        area=area,
        d=d,
        d2=d2,
        e=e,
        f=f,
        a=a,
        b=b,
        storage_coef=s_coef,
        alpha=alpha,
        storage=[store],
        quickflow=[quick],
        slowflow=[slow]
        )
end


function update_state(s_node::IHACRESNode, storage::Float64, e_rainfall::Float64, et::Float64, qflow_store::Float64, sflow_store::Float64, outflow::Float64, level::Float64)
    push!(s_node.storage, storage)
    push!(s_node.effective_rainfall, e_rainfall)
    push!(s_node.et, et)
    push!(s_node.outflow, outflow)

    push!(s_node.quickflow, qflow_store)
    push!(s_node.slowflow, sflow_store)
    push!(s_node.level, level)
end


"""
Run node to calculate outflow and update state.

Parameters
----------
timestep: int, time step
rain: float, rainfall
evap: float, evapotranspiration
extractions: float, irrigation and other water extractions
gw_exchange: float, flux in ML - positive is contribution to stream, negative is infiltration
loss: float,

Returns
----------
float, outflow from node
"""
function run_node!(s_node::IHACRESNode,
                  rain::Float64,
                  evap::Float64,
                  inflow::Float64,
                  ext::Float64,
                  gw_exchange::Float64=0.0,
                  loss::Float64=0.0)::Tuple{Float64, Float64}
    arr_len = length(s_node.storage)
    current_store = s_node.storage[arr_len]

    interim_results = [0.0, 0.0, 0.0]
    @ccall ihacres.calc_ft_interim(interim_results::Ptr{Cdouble},
                                   current_store::Cdouble,
                                   rain::Cdouble,
                                   s_node.d::Cdouble,
                                   s_node.d2::Cdouble,
                                   s_node.alpha::Cdouble)::Cvoid

    (mf, e_rainfall, recharge) = interim_results

    et::Float64 = @ccall ihacres.calc_ET(
        s_node.e::Cdouble,
        evap::Cdouble,
        mf::Cdouble,
        s_node.f::Cdouble,
        s_node.d::Cdouble
    )::Cdouble

    cmd::Float64 = @ccall ihacres.calc_cmd(
        current_store::Cdouble,
        rain::Cdouble,
        s_node.d::Cdouble,
        s_node.d2::Cdouble,
        s_node.alpha::Cdouble
    )::Cdouble

    # var inflow = 0.0
    # for nid in s_node.prev_node:
    #     inflow += s_node.prev_node[nid].run(timestep, rain_evap, ext)
    # # End for
    push!(s_node.inflow, inflow)

    flow_results = [0.0, 0.0, 0.0]
    @ccall ihacres.calc_ft_flows(
        flow_results::Ptr{Cdouble},
        s_node.quickflow[arr_len]::Cdouble,
        s_node.slowflow[arr_len]::Cdouble,
        e_rainfall::Cdouble,
        recharge::Cdouble,
        s_node.area::Cdouble,
        s_node.a::Cdouble,
        s_node.b::Cdouble,
        loss::Cdouble
    )::Cvoid

    (quick_store, slow_store, outflow) = flow_results

    # if self.next_node:  # and ('dam' not in self.next_node.node_type):
    #     cmd, outflow = routing(cmd, s_node.storage_coef, inflow, outflow, ext, gamma=gw_exchange)
    # else:
    #     outflow = calc_outflow(outflow, ext)
    # # End if
    routing_res = [0.0, 0.0]
    @ccall ihacres.routing(
        routing_res::Ptr{Cdouble},
        cmd::Cdouble,
        s_node.storage_coef::Cdouble,
        inflow::Cdouble,
        outflow::Cdouble,
        ext::Cdouble,
        gw_exchange::Cdouble)::Cvoid

    (cmd, outflow) = routing_res

    level::Float64 = @ccall ihacres.calc_ft_level(
        outflow::Cdouble,
        s_node.level_params::Ptr{Cdouble}
    )::Cdouble

    update_state(s_node, cmd, e_rainfall, et, quick_store, slow_store, outflow, level)

    return outflow, level
end


"""
"""
function update_params!(node::IHACRESNode, area::Float64, d::Float64, d2::Float64, e::Float64, f::Float64,
                        a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                        p1::Float64, p2::Float64, p3::Float64, p4::Float64, p5::Float64, p6::Float64, p7::Float64, p8::Float64, CTF::Float64)
    node.area = area
    node.d = d
    node.d2 = d2
    node.e = e
    node.f = f
    node.a = a
    node.b = b
    node.storage_coef = s_coef
    node.alpha = alpha
    node.level_params = [p1, p2, p3, p4, p5, p6, p7, p8, CTF]

    return nothing
end



function reset!(s_node::IHACRESNode)
    s_node.storage = [s_node.storage[1]]

    s_node.quickflow = [s_node.quickflow[1]]
    s_node.slowflow = [s_node.slowflow[1]]
    s_node.outflow = []
    s_node.effective_rainfall = []
    s_node.et = []
    s_node.inflow = []
end