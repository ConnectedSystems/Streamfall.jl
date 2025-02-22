"""
Portions of the GR4J implementation is adapted from Python code written by Andrew MacDonald (2014).

Per license requirements, the full license conditions included below.


Copyright (c) 2014, Andrew MacDonald
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.

    * The name of the author may not be used to
    endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
"""

using Parameters
using ModelParameters


abstract type GRNJNode <: NetworkNode end


"""
GR4J Node

A four-parameter model with two stores.

# Parameters
- `x1` : maximum capacity of the production store (mm) (> 0)
- `x2` : groundwater exchange coefficient (mm) (value < and > 0 possible)
- `x3` : one day ahead maximum capacity of the routing store (mm, > 0)
- `x4` : time base of unit hydrograph UH1 (days, > 0.5)

# References
1. Perrin, C., Michel, C., Andréassian, V., 2003.
    Improvement of a parsimonious model for streamflow simulation.
    Journal of Hydrology 279, 275-289.
    https://doi.org/10.1016/S0022-1694(03)00225-7

2. MacDonald, A. 2014.
    Python GR4J
    https://github.com/amacd31/gr4j
"""
Base.@kwdef mutable struct GR4JNode{P,A<:AbstractFloat} <: GRNJNode
    const name::String
    const area::A

    # Parameters
    # x1 : maximum capacity of the production store (mm) (> 0)
    # x2 : groundwater exchange coefficient (mm) (value < and > 0 possible)
    # x3 : one day ahead maximum capacity of the routing store (mm, > 0)
    # x4 : time base of unit hydrograph UH1 (days, > 0.5)
    X1::P = Param(350.0, bounds=(1.0, 1500.0))
    X2::P = Param(0.0, bounds=(-10.0, 5.0))
    X3::P = Param(40.0, bounds=(1.0, 500.0))
    X4::P = Param(0.5, bounds=(0.5, 10.0))

    # stores
    p_store::Vector{A} = [0.0]
    r_store::Vector{A} = [0.0]

    UH1::Vector{Vector{A}} = []
    UH2::Vector{Vector{A}} = []

    uh1_ordinates::Vector{A} = []
    uh2_ordinates::Vector{A} = []

    outflow::Vector{A} = []

    obj_func::Function = obj_func
end


function prep_state!(node::GR4JNode, timesteps::Int64)
    resize!(node.p_store, timesteps + 1)
    resize!(node.r_store, timesteps + 1)
    node.p_store[2:end] .= 0.0
    node.r_store[2:end] .= 0.0

    # Prep cache
    X4 = node.X4.val
    nUH1 = Int(ceil(X4))
    nUH2 = Int(ceil(2.0 * X4))
    cUH1, cUH2 = fill(0.0, nUH1), fill(0.0, nUH2)
    node.UH1 = fill(cUH1, timesteps)
    node.UH2 = fill(cUH2, timesteps)
    node.uh1_ordinates = Vector{Float64}(undef, nUH1)
    node.uh2_ordinates = Vector{Float64}(undef, nUH2)

    node.outflow = fill(0.0, timesteps)
end

function GR4JNode(name::String, spec::AbstractDict)
    n = create_node(GR4JNode, name, spec["area"])
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
    run_timestep!(
        node::GR4JNode, climate::Climate, timestep::Int;
        inflow=nothing, extraction=nothing, exchange=nothing
    )
    run_timestep!(
        node::GR4JNode, rain::Float64, et::Float64, ts::Int64;
        inflow=nothing, extraction=nothing, exchange=nothing
    )

Run given GR4J node for a time step.
"""
function run_timestep!(
    node::GR4JNode, climate::Climate, timestep::Int;
    inflow=nothing, extraction=nothing, exchange=nothing
)
    P, E = climate_values(node, climate, timestep)

    run_timestep!(node, P, E, timestep; inflow=inflow, extraction=extraction, exchange=exchange)
end
function run_timestep!(
    node::GR4JNode, rain::Float64, et::Float64, ts::Int64;
    inflow=nothing, extraction=nothing, exchange=nothing
)
    uh1_cache = node.uh1_ordinates
    uh2_cache = node.uh2_ordinates
    res = run_gr4j(
        rain, et, node.X1.val, node.X2.val, node.X3.val, node.X4.val, node.area,
        node.UH1[ts], node.UH2[ts], uh1_cache, uh2_cache;
        p_store=node.p_store[ts], r_store=node.r_store[ts]
    )
    Q, p_s, r_s, UH1, UH2 = res

    node_name = node.name
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)

    if !isnothing(inflow)
        Q = Q + in_flow + ex - wo
    end

    update_state!(node, ts, p_s, r_s, Q, UH1, UH2)

    return Q
end

"""
    update_state!(node::GR4JNode, ps, rs, q, UH1, UH2)::Nothing
    update_state!(node::GR4JNode, ts::Int64, ps, rs, q, UH1, UH2)::Nothing

Update GR4J node state.
"""
function update_state!(node::GR4JNode, ps, rs, q, UH1, UH2)::Nothing
    append!(node.p_store, ps)
    append!(node.r_store, rs)
    append!(node.outflow, q)
    push!(node.UH1, UH1)
    push!(node.UH2, UH2)

    return nothing
end
function update_state!(node::GR4JNode, ts::Int64, ps, rs, q, UH1, UH2)::Nothing
    node.p_store[ts+1] = ps
    node.r_store[ts+1] = rs
    node.outflow[ts] = q

    node.UH1[ts] = UH1
    node.UH2[ts] = UH2

    return nothing
end

"""
    update_params!(node::GR4JNode, X1::Float64, X2::Float64, X3::Float64, X4::Float64)::Nothing

Update parameters for GR4J.
"""
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
    # Stores
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
function s_curve(t::F, x4::F; uh2::Bool=false)::F where {F<:Float64}
    if t <= 0.0
        return 0.0
    end

    ordinate::F = 0.0
    if t < x4
        ordinate = (t / x4)^2.5
    else
        if uh2 && (t < 2 * x4)
            ordinate = 1.0 - 0.5 * (2 - t / x4)^2.5
        else
            # t >= x4 if uh1, or
            # t >= 2*x4 if uh2
            ordinate = 1.0
        end
    end

    return ordinate
end


"""
    run_gr4j(
        P::F, E::F, X1::F, X2::F, X3::F, X4::F, area::F,
        UH1::Vector{Float64}, UH2::Vector{Float64},
        uh1_ordinates::Vector{Float64}, uh2_ordinates::Vector{Float64};
        p_store=0.0, r_store=0.0
    )::Tuple where {F<:Float64}

Generated simulated streamflow for given rainfall and potential evaporation.

# Parameters
- `P` : Catchment average rainfall
- `E` : Catchment average potential evapotranspiration
- `X1` : Maximum capacity of production store (in mm; > 0)
- `X2` : Groundwater exchange coefficient (in mm; value < and > 0 possible)
- `X3` : Maximum capacity of routing store (in mm; > 0)
- `X4` : Time base of the unit hydrograph (in days, > 0.5)
- `area` : Catchment area
- `UH1` : Quickflow store
- `UH2` : Baseflow store
- `uh1_ordinates` : The proportion of rainfall converted to quickflow for each timestep
- `uh2_ordinates` : The proportion of rainfall converted to slowflow for each timestep
- `p_store` : Initial production store
- `r_store` : Initial state store

# Returns
Tuple of:
- Simulated outflow [ML/day]
- intermediate states:
  - p_store (initial production / percolation)
  - r_store (initial state)
  - UH1 (Quickflow)
  - UH2 (Slowflow)
"""
function run_gr4j(
    P::F, E::F, X1::F, X2::F, X3::F, X4::F, area::F,
    UH1::Vector{F}, UH2::Vector{F},
    uh1_ordinates::Vector{F}, uh2_ordinates::Vector{F};
    p_store=0.0, r_store=0.0
)::Tuple where {F<:Float64}
    nUH1::Int64 = Int(ceil(X4))
    nUH2::Int64 = Int(ceil(2.0 * X4))

    @inbounds for t in 2:(nUH1+1)
        t_f = Float64(t)
        uh1_ordinates[t-1] = s_curve(t_f, X4) - s_curve(t_f - 1.0, X4)
    end

    @inbounds for t in 2:(nUH2+1)
        t_f = Float64(t)
        uh2_ordinates[t-1] = s_curve(t_f, X4, uh2=true) - s_curve(t_f - 1.0, X4, uh2=true)
    end

    Q::F = 0.0
    if P > E
        net_evap = 0.0
        scaled_net_precip = min((P - E) / X1, 13.0)
        tanh_scaled_net_precip = tanh(scaled_net_precip)

        tmp_a = (p_store / X1)^2
        numer = (X1 * (1.0 - tmp_a) * tanh_scaled_net_precip)
        denom = (1.0 + p_store / X1 * tanh_scaled_net_precip)
        reservoir_production = numer / denom

        routed_volume = P - E - reservoir_production
    else
        scaled_net_evap = min((E - P) / X1, 13.0)
        tanh_scaled_net_evap = tanh(scaled_net_evap)

        ps_div_x1 = (2.0 - p_store / X1) * tanh_scaled_net_evap
        net_evap = (p_store * (ps_div_x1) /
                    (1.0 + (1.0 - p_store / X1) * tanh_scaled_net_evap))

        reservoir_production = 0.0
        routed_volume = 0.0
    end

    p_store = p_store - net_evap + reservoir_production

    tmp_a::F = (p_store / 2.25 / X1)^4
    tmp_b::F = (1 + tmp_a)^0.25
    percolation = p_store / tmp_b

    routed_volume = routed_volume + (p_store - percolation)

    @inbounds for i in 1:nUH1-1
        UH1[i] = UH1[i+1] + uh1_ordinates[i] * routed_volume
    end
    UH1[end] = uh1_ordinates[end] * routed_volume

    @inbounds for j in 1:nUH2-1
        UH2[j] = UH2[j+1] + uh2_ordinates[j] * routed_volume
    end
    UH2[end] = uh2_ordinates[end] * routed_volume

    tmp_a = (r_store / X3)^3.5
    groundwater_exchange::F = X2 * tmp_a
    r_store = max(0.0, r_store + UH1[1] * 0.9 + groundwater_exchange)

    tmp_a = (r_store / X3)^4
    tmp_b = (1 + tmp_a)^0.25
    R2 = r_store / tmp_b
    QR = r_store - R2
    r_store = R2
    QD = max(0.0, UH2[1] * 0.1 + groundwater_exchange)

    Q = (QR + QD) * area

    return Q, percolation, r_store, UH1, UH2
end
