using Parameters
using ModelParameters


abstract type GR4JNode <: NetworkNode end


Base.@kwdef mutable struct SimpleGR4JNode{A <: Union{Param, Real}} <: GR4JNode
    @network_node

    # parameters
    X1::A = Param(303.6276, bounds=(50.0, 500.0))
    X2::A = Param(0.3223, bounds=(0.0, 1.0))
    X3::A = Param(6.4975, bounds=(0.1, 10.0))
    X4::A = Param(0.2948, bounds=(0.1, 0.999))

    # stores
    p_store::Float64 = 0.0
    r_store::Float64 = 0.0

    UH1::Array{Array{Float64}} = []
    UH2::Array{Array{Float64}} = []

    outflow::Array{Float64} = []
end


function run_node!(node::GR4JNode, climate::Climate)
    timesteps = sim_length(climate)
    for ts in 1:timesteps
        run_node!(node, climate, ts)
    end

    return node.outflow
end


function run_node!(node::GR4JNode, climate::Climate, timestep::Int)
    P, E = climate_values(node, climate, timestep)

    res = run_gr4j(P, E, node.X1, node.X2, node.X3, node.X4, node.p_store, node.r_store)
    Q, p_s, r_s, UH1, UH2 = res

    update_state!(node, p_s, r_s, Q, UH1, UH2)
end


function update_state!(node::GR4JNode, ps, rs, q, UH1, UH2)
    node.p_store = ps
    node.r_store = rs
    append!(node.outflow, q)
    append!(node.UH1, UH1)
    append!(node.UH2, UH2)
end


"""
    s_curve(t::Float64, x4::Float64, uh2::Bool = false)::Float64

Determine unit hydrograph ordinates.
"""
function s_curve(t::Float64, x4::Float64; uh2::Bool = false)::Float64
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
Generated simulated streamflow for given rainfall and potential evaporation.

# Parameters
- P: Catchment average rainfall
- E: Catchment average potential evapotranspiration
- X1 - X4 : X parameters for the model
- p_store : initial production store
- r_store : initial state store

# Returns
- tuple of simulated outflow, and intermediate states: p_store, r_store, UH1, UH2
"""
function run_gr4j(P::Float64, E::Float64,
                  X1::Float64, X2::Float64, X3::Float64, X4::Float64,
                  p_store::Float64=0.0, r_store::Float64=0.0)::Tuple
    nUH1 = int(ceil(X4))
    nUH2 = int(ceil(2.0*X4))

    uh1_ordinates = Array{Float64}(undef, nUH1)
    uh2_ordinates = Array{Float64}(undef, nUH2)

    UH1, UH2 = Array{Float64}(0.0, nUH1), Array{Float64}(0.0, nUH2)

    for t in 2:nUH1
        t_f = float(t)
        uh1_ordinates[t - 1] = s_curve(t_f, X4) - s_curve(t_f, X4)
        uh2_ordinates[t - 1] = s_curve(t_f, X4, uh2=true) - s_curve(t_f-1, X4, uh2=true)
    end
    
    for t in (nUH1+2):nUH2
        t_f = float(t)
        uh2_ordinates[t - 1] = s_curve(t_f, X4, uh2=true) - s_curve(t_f-1, X4, uh2=true)
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

    p_store = p_store - net_evap + reservoir_production

    tmp_a = (p_store/2.25/X1)^4
    tmp_b = (1 + tmp_a)^0.25
    percolation = p_store / tmp_b

    routed_volume = routed_volume + (p_store-percolation)
    p_store = percolation

    # Check these loops too
    for i in 1:(length(UH1)-1)
        UH1[i] = UH1[i+1] + uh1_ordinates[i]*routing_pattern
    end
    UH1[end] = uh1_ordinates[end] * routing_pattern

    for j in 1:(length(UH2)-2)
        UH2[j] = UH2[j+1] + uh2_ordinates[j]*routing_pattern
    end
    UH2[end] = uh2_ordinates[end] * routing_pattern

    tmp_a = (r_store / X3)^3.5
    groundwater_exchange = X2 * tmp_a
    r_store = max(0.0, r_store + UH1[1] * 0.9 + groundwater_exchange)

    tmp_a = (r_store / X3)^4
    tmp_b = (1 + tmp_a)^0.25
    R2 = r_store / tmp_b
    QR = r_store - R2
    r_store = R2
    QD = max(0, UH2[1]*0.1+groundwater_exchange)
    Q = QR + QD

    return Q, p_store, r_store, UH1, UH2
end
