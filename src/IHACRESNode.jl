using Parameters
using ModelParameters


abstract type IHACRESNode <: NetworkNode end


Base.@kwdef mutable struct BilinearNode{A <: Union{Param, Real}} <: IHACRESNode
    @network_node

    # https://wiki.ewater.org.au/display/SD41/IHACRES-CMD+-+SRG
    d::A = Param(200.0, bounds=(10.0, 550.0))  # flow threshold
    d2::A = Param(2.0, bounds=(0.0001, 10.0))   # flow threshold2
    e::A = Param(1.0, bounds=(0.1, 1.5))  # temperature to PET conversion factor
    f::A = Param(0.8, bounds=(0.01, 3.0))  # plant stress threshold factor (multiplicative factor of d)
    a::A = Param(0.9, bounds=(0.1, 100.0))
    b::A = Param(0.1, bounds=(0.0, 1.0))
    storage_coef::A = Param(2.9, bounds=(0.2, 10.0))
    alpha::A = Param(0.1, bounds=(1e-5, 1 - 1/10^9))

    level_params::Array{A, 1} = [
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

    storage::Array{Float64} = [100.0]
    quick_store::Array{Float64} = [0.0]
    slow_store::Array{Float64} = [0.0]
    outflow::Array{Float64} = []
    effective_rainfall::Array{Float64} = []
    et::Array{Float64} = []
    inflow::Array{Float64} = []
    level::Array{Float64} = []
end


function BilinearNode(node_id::String, spec::Dict)
    route = true
    if isnothing(spec["inlets"]) || isempty(spec["inlets"])
        route = false
    end
        
    n = BilinearNode{Param}(; node_id=node_id, area=spec["area"], 
                           route=route)

    node_params = spec["parameters"]
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


function BilinearNode(node_id::String, area::Float64, route::Bool, d::Float64, d2::Float64, e::Float64, f::Float64, 
                    a::Float64, b::Float64, s_coef::Float64, alpha::Float64, 
                    store::Float64, quick::Float64, slow::Float64)
    return BilinearNode{Float64}(
        node_id=node_id,
        area=area,
        route=route,
        d=d,
        d2=d2,
        e=e,
        f=f,
        a=a,
        b=b,
        storage_coef=s_coef,
        alpha=alpha,
        storage=[store],
        quick_store=[quick],
        slow_store=[slow]
    )
end


function update_state(s_node::IHACRESNode, storage::Float64, e_rainfall::Float64, et::Float64, qflow_store::Float64, sflow_store::Float64, outflow::Float64, level::Float64)
    push!(s_node.storage, storage)
    push!(s_node.effective_rainfall, e_rainfall)
    push!(s_node.et, et)
    push!(s_node.outflow, outflow)

    push!(s_node.quick_store, qflow_store)
    push!(s_node.slow_store, sflow_store)
    push!(s_node.level, level)
end


"""
Run node with ET data to calculate outflow and update state.

Parameters
----------
s_node
rain: float, rainfall
evap: float, evapotranspiration
inflow: float, inflow from previous node
ext: float, irrigation and other water extractions
gw_exchange: float, flux in ML - positive is contribution to stream, negative is infiltration
loss: float,
current_store: replacement cmd value
quick_store: replacement quick_store value
slow_store: replacement slow_store value

Returns
----------
float, outflow from node
"""
function run_node!(s_node::BilinearNode,
                   rain::Float64,
                   evap::Float64,
                   inflow::Float64,
                   ext::Float64,
                   gw_exchange::Float64=0.0,
                   loss::Float64=0.0;
                   current_store=nothing,
                   quick_store=nothing,
                   slow_store=nothing)::Tuple{Float64, Float64}
    if !isnothing(current_store)
        s_node.storage[end] = current_store
        # current_store = s_node.storage[end]
    end

    if !isnothing(quick_store)
        s_node.quick_store[end] = quick_store
        # quick_store = s_node.quick_store[end]
    end

    if !isnothing(slow_store)
        s_node.slow_store[end] = slow_store
        # slow_store = s_node.slow_store[end]
    end

    current_store = s_node.storage[end]
    quick_store = s_node.quick_store[end]
    slow_store = s_node.slow_store[end]

    interim_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_interim(interim_results::Ptr{Cdouble},
                                   current_store::Cdouble,
                                   rain::Cdouble,
                                   s_node.d::Cdouble,
                                   s_node.d2::Cdouble,
                                   s_node.alpha::Cdouble)::Cvoid

    # @assert any(isnan.(interim_results)) == false
    (mf, e_rainfall, recharge) = interim_results

    et::Float64 = @ccall IHACRES.calc_ET(
        s_node.e::Cdouble,
        evap::Cdouble,
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
        s_node.b::Cdouble,
        loss::Cdouble
    )::Cvoid

    # @assert any(isnan.(flow_results)) == false
    (nq_store, ns_store, outflow) = flow_results

    # if self.next_node:  # and ('dam' not in self.next_node.node_type):
    #     cmd, outflow = routing(cmd, s_node.storage_coef, inflow, outflow, ext, gamma=gw_exchange)
    # else:
    #     outflow = calc_outflow(outflow, ext)
    # # End if
    if s_node.route
        routing_res = [0.0, 0.0]
        @ccall IHACRES.routing(
            routing_res::Ptr{Cdouble},
            cmd::Cdouble,
            s_node.storage_coef::Cdouble,
            inflow::Cdouble,
            outflow::Cdouble,
            ext::Cdouble,
            gw_exchange::Cdouble)::Cvoid

        # @assert any(isnan.(routing_res)) == false
        (vol, outflow) = routing_res
    end

    level_params = Array{Float64}(s_node.level_params)
    level::Float64 = @ccall IHACRES.calc_ft_level(
        outflow::Cdouble,
        level_params::Ptr{Cdouble}
    )::Cdouble

    push!(s_node.inflow, inflow)
    update_state(s_node, cmd, e_rainfall, et, nq_store, ns_store, outflow, level)

    return outflow, level
end


"""
Extract node parameter values and bounds
"""
function param_info(node::IHACRESNode; with_level::Bool = true)::Tuple
    tmp = Model(node)
    values = collect(tmp.val)
    bounds = collect(tmp.bounds)

    if with_level
        level_param_vals = map(x->x.val, node.level_params)
        append!(values, level_param_vals)

        level_param_bounds = map(x->x.bounds, node.level_params)
        append!(bounds, level_param_bounds)
    end
    
    return values, bounds
end


"""
Run node with temperature data to calculate outflow and update state.

Parameters
----------
s_node
rain: float, rainfall
temp: float, temperature
inflow: float, inflow from previous node
ext: float, irrigation and other water extractions
gw_exchange: float, flux in ML - positive is contribution to stream, negative is infiltration
loss: float,
current_store: replacement cmd value
quick_store: replacement quick_store value
slow_store: replacement slow_store value

Returns
----------
float, outflow from node
"""
function run_node_with_temp!(s_node::BilinearNode,
                   rain::Float64,
                   temp::Float64,
                   inflow::Float64,
                   ext::Float64,
                   gw_exchange::Float64=0.0,
                   loss::Float64=0.0;
                   current_store=nothing,
                   quick_store=nothing,
                   slow_store=nothing)::Tuple{Float64, Float64}
    if !isnothing(current_store)
        s_node.storage[end] = current_store
        # current_store = s_node.storage[end]
    end

    if !isnothing(quick_store)
        s_node.quick_store[end] = quick_store
        # quick_store = s_node.quick_store[end]
    end

    if !isnothing(slow_store)
        s_node.slow_store[end] = slow_store
        # slow_store = s_node.slow_store[end]
    end

    current_store = s_node.storage[end]
    quick_store = s_node.quick_store[end]
    slow_store = s_node.slow_store[end]

    interim_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_interim(interim_results::Ptr{Cdouble},
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

    # convert to areal average
    # et = et / s_node.area

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
        s_node.b::Cdouble,
        loss::Cdouble
    )::Cvoid

    @assert any(isnan.(flow_results)) == false
    (nq_store, ns_store, outflow) = flow_results

    # if self.next_node:  # and ('dam' not in self.next_node.node_type):
    #     cmd, outflow = routing(cmd, s_node.storage_coef, inflow, outflow, ext, gamma=gw_exchange)
    # else:
    #     outflow = calc_outflow(outflow, ext)
    # # End if
    if s_node.route
        routing_res = [0.0, 0.0]
        @ccall IHACRES.routing(
            routing_res::Ptr{Cdouble},
            cmd::Cdouble,
            s_node.storage_coef::Cdouble,
            inflow::Cdouble,
            outflow::Cdouble,
            ext::Cdouble,
            gw_exchange::Cdouble)::Cvoid

        @assert any(isnan.(routing_res)) == false
        (vol, outflow) = routing_res
    end

    level_params = Array{Float64}(s_node.level_params)
    level::Float64 = @ccall IHACRES.calc_ft_level(
        outflow::Cdouble,
        level_params::Ptr{Cdouble}
    )::Cdouble

    push!(s_node.inflow, inflow)
    update_state(s_node, cmd, e_rainfall, et, nq_store, ns_store, outflow, level)

    return outflow, level
end


"""
"""
function update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,
                        a::Float64, b::Float64, s_coef::Float64, alpha::Float64)::Nothing
    node.d = Param(d, bounds=node.d.bounds)
    node.d2 = Param(d2, bounds=node.d2.bounds)
    node.e = Param(e, bounds=node.e.bounds)
    node.f = Param(f, bounds=node.f.bounds)
    node.a = Param(a, bounds=node.a.bounds)
    node.b = Param(b, bounds=node.b.bounds)
    node.storage_coef = Param(s_coef, bounds=node.storage_coef.bounds)
    node.alpha = Param(alpha, bounds=node.alpha.bounds)

    return nothing
end


"""
Update all parameters
"""
function update_params!(node::BilinearNode{Param}, d::Float64, d2::Float64, e::Float64, f::Float64,
                        a::Float64, b::Float64, s_coef::Float64, alpha::Float64,
                        p1::Float64, p2::Float64, p3::Float64, p4::Float64, p5::Float64, p6::Float64, p7::Float64, p8::Float64, CTF::Float64)::Nothing
    node.d = Param(d, bounds=node.d.bounds)
    node.d2 = Param(d2, bounds=node.d2.bounds)
    node.e = Param(e, bounds=node.e.bounds)
    node.f = Param(f, bounds=node.f.bounds)
    node.a = Param(a, bounds=node.a.bounds)
    node.b = Param(b, bounds=node.b.bounds)
    node.storage_coef = Param(s_coef, bounds=node.storage_coef.bounds)
    node.alpha = Param(alpha, bounds=node.alpha.bounds)

    n_lparams = node.level_params
    node.level_params = [
        Param(p1, bounds=n_lparams[1].bounds)
        Param(p2, bounds=n_lparams[2].bounds)
        Param(p3, bounds=n_lparams[3].bounds)
        Param(p4, bounds=n_lparams[4].bounds)
        Param(p5, bounds=n_lparams[5].bounds)
        Param(p6, bounds=n_lparams[6].bounds)
        Param(p7, bounds=n_lparams[7].bounds)
        Param(p8, bounds=n_lparams[8].bounds)
        Param(CTF, bounds=n_lparams[9].bounds)
    ] 

    return nothing
end



function reset!(s_node::IHACRESNode)::Nothing
    s_node.storage = [s_node.storage[1]]
    s_node.quick_store = [s_node.quick_store[1]]
    s_node.slow_store = [s_node.slow_store[1]]
    s_node.outflow = []
    s_node.level = []
    s_node.effective_rainfall = []
    s_node.et = []
    s_node.inflow = []

    return nothing
end