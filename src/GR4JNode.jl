using Parameters
using ModelParameters


abstract type GRNJNode <: NetworkNode end


Base.@kwdef mutable struct GR4JNode{A <: Union{Param, Real}} <: GRNJNode
    @network_node

    # parameters
    X1::A = Param(350.0, bounds=(1.0, 1500.0))
    X2::A = Param(0.0, bounds=(-10.0, 5.0))
    X3::A = Param(40.0, bounds=(1.0, 500.0))
    X4::A = Param(0.5, bounds=(0.5, 10.0))

    # x1 : maximum capacity of the production store (mm) (> 0)
    # x2 : groundwater exchange coefficient (mm) (value < and > 0 possible)
    # x3 : one day ahead maximum capacity of the routing store (mm, > 0)
    # x4 : time base of unit hydrograph UH1 (days, > 0.5)

    # stores
    p_store::Array{Float64} = [0.0]
    r_store::Array{Float64} = [0.0]

    UH1::Array{Array{Float64}} = []
    UH2::Array{Array{Float64}} = []

    outflow::Array{Float64} = []
end


# function prep_state!(node::HyModNode, sim_length::Int64)
#     node.UH1 = fill(undef, sim_length)
#     node.UH2 = fill(undef, sim_length)
#     node.outflow = fill(undef, sim_length)
# end


function GR4JNode(name::String, spec::Dict)
    n = GR4JNode{Param}(; name=name, area=spec["area"])
    node_params = spec["parameters"]
    node_params["p_store"] = [node_params["initial_p_store"]]
    node_params["r_store"] = [node_params["initial_r_store"]]

    delete!(node_params, "initial_p_store")
    delete!(node_params, "initial_r_store")

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
    run_node!(node::GR4JNode, climate::Climate)

Run GR4J node for all time steps in given climate sequence.
"""
function run_node!(node::GR4JNode, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    # prep_state!(node, sim_length)
    for ts in 1:timesteps
        run_node!(node, climate, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow
end


"""
    run_node!(node::GR4JNode, climate::Climate, timestep::Int;
              inflow::Float64, extraction::Float64, exchange::Float64)

Run given GR4J node for a time step.
"""
function run_node!(node::GR4JNode, climate::Climate, timestep::Int; inflow=nothing, extraction=nothing, exchange=nothing)
    P, E = climate_values(node, climate, timestep)

    res = run_gr4j(P, E, node.X1, node.X2, node.X3, node.X4, node.area, node.p_store[end], node.r_store[end])
    Q, p_s, r_s, UH1, UH2 = res

    ts = timestep
    node_name = node.name
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)

    if !isnothing(inflow)
        Q = Q + in_flow + ex - wo
    end

    update_state!(node, p_s, r_s, Q, UH1, UH2)
end


function update_state!(node::GR4JNode, ps, rs, q, UH1, UH2)
    append!(node.p_store, ps)
    append!(node.r_store, rs)
    append!(node.outflow, q)
    push!(node.UH1, UH1)
    push!(node.UH2, UH2)
end


function update_params!(node::GR4JNode, X1::Float64, X2::Float64, X3::Float64, X4::Float64)::Nothing
    node.X1 = Param(X1, bounds=node.X1.bounds)
    node.X2 = Param(X2, bounds=node.X2.bounds)
    node.X3 = Param(X3, bounds=node.X3.bounds)
    node.X4 = Param(X4, bounds=node.X4.bounds)

    return nothing
end


"""
    reset!(node::GR4JNode)::Nothing

Reset node. Clears all states back to their initial values.
"""
function reset!(node::GR4JNode)::Nothing
    # stores
    node.p_store = [node.p_store[1]]
    node.r_store = [node.r_store[1]]

    node.UH1 = []
    node.UH2 = []

    node.outflow = []

    return nothing
end


"""
    s_curve(t::Float64, x4::Float64, uh2::Bool = false)::Float64

Determine unit hydrograph ordinates.
"""
function s_curve(t, x4; uh2::Bool = false)::Float64
    if t <= 0.0
        return 0.0
    end

    ordinate::Float64 = 0.0
    if t < x4
        ordinate = (t/x4)^2.5
    else
        if uh2 && (t < 2*x4)
            ordinate = 1.0 - 0.5*(2 - t/x4)^2.5
        else
            # t >= x4 if uh1, or
            # t >= 2*x4 if uh2
            ordinate = 1.0
        end
    end

    return ordinate
end


"""
    run_gr4j(P::Float64, E::Float64,
             X1::Float64, X2::Float64, X3::Float64, X4::Float64, area::Float64,
             p_store::Float64=0.0, r_store::Float64=0.0)::Tuple

Generated simulated streamflow for given rainfall and potential evaporation.

# Parameters
- P : Catchment average rainfall
- E : Catchment average potential evapotranspiration
- X1 - X4 : X parameters
- area : Catchment area
- p_store : Initial production store
- r_store : Initial state store

# Returns
- tuple of simulated outflow [ML/day], and intermediate states: p_store, r_store, UH1, UH2
"""
function run_gr4j(P, E, X1, X2, X3, X4, area, p_store=0.0, r_store=0.0)::Tuple
    nUH1 = Int(ceil(X4))
    nUH2 = Int(ceil(2.0*X4))

    uh1_ordinates = Array{Float64}(undef, nUH1)
    uh2_ordinates = Array{Float64}(undef, nUH2)

    UH1, UH2 = fill(0.0, nUH1), fill(0.0, nUH2)

    for t in 2:(nUH1+1)
        t_f = Float64(t)
        uh1_ordinates[t - 1] = s_curve(t_f, X4) - s_curve(t_f-1.0, X4)
    end

    for t in 2:(nUH2+1)
        t_f = Float64(t)
        uh2_ordinates[t - 1] = s_curve(t_f, X4, uh2=true) - s_curve(t_f-1.0, X4, uh2=true)
    end

    Q = 0.0
    if P > E
        net_evap = 0.0
        scaled_net_precip = min((P - E)/X1, 13.0)
        tanh_scaled_net_precip = tanh(scaled_net_precip)

        tmp_a = (p_store/X1)^2
        numer = (X1 * (1.0 - tmp_a) * tanh_scaled_net_precip)
        denom = (1.0 + p_store/X1 * tanh_scaled_net_precip)
        reservoir_production = numer / denom

        routed_volume = P-E-reservoir_production
    else
        scaled_net_evap = min((E - P)/X1, 13.0)
        tanh_scaled_net_evap = tanh(scaled_net_evap)

        ps_div_x1 = (2.0 - p_store/X1) * tanh_scaled_net_evap
        net_evap = (p_store * (ps_div_x1) /
                    (1.0 + (1.0 - p_store/X1) * tanh_scaled_net_evap))

        reservoir_production = 0.0
        routed_volume = 0.0
    end

    p_store = p_store - net_evap + reservoir_production

    tmp_a = (p_store / 2.25 / X1)^4
    tmp_b = (1 + tmp_a)^0.25
    percolation = p_store / tmp_b

    routed_volume = routed_volume + (p_store-percolation)
    p_store = percolation

    # Check these loops too
    for i in 1:(length(UH1)-1)
        UH1[i] = UH1[i+1] + uh1_ordinates[i]*routed_volume
    end
    UH1[end] = uh1_ordinates[end] * routed_volume

    for j in 1:(length(UH2)-1)
        UH2[j] = UH2[j+1] + uh2_ordinates[j]*routed_volume
    end
    UH2[end] = uh2_ordinates[end] * routed_volume

    tmp_a = (r_store / X3)^3.5
    groundwater_exchange = X2 * tmp_a
    r_store = max(0.0, r_store + UH1[1] * 0.9 + groundwater_exchange)

    tmp_a = (r_store / X3)^4
    tmp_b = (1 + tmp_a)^0.25
    R2 = r_store / tmp_b
    QR = r_store - R2
    r_store = R2
    QD = max(0.0, UH2[1]*0.1+groundwater_exchange)
    Q = (QR + QD) * area * 1000.0 / 86.4
    # (QR + QD) * area * 1000.0/86400.0  # m³/sec

    return Q, p_store, r_store, UH1, UH2
end
