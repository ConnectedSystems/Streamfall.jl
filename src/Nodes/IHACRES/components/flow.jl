approx_zero(x) = x + one(x) ≈ one(x)

"""
    convert_taus(tau::Float64, v::Float64)

Convert τ to α and β

# Arguments
- `tau` : time constant τ
- `v`   : proportion of flow v (between 0 and 1)

# Returns
Tuple of (alpha, beta)
"""
@inline function convert_taus(tau::Float64, v::Float64)::Tuple{Float64,Float64}
    return (exp(-1.0 / tau), v * (1.0 - alpha))
end

"""
    calc_stores(tau::Float64, prev_store::Float64, v::Float64, e_rainfall::Float64)

Calculate flow stores.
"""
function calc_stores(tau::Float64, prev_store::Float64, v::Float64, e_rainfall::Float64)::Float64
    alpha, beta = convert_taus(tau, v)
    store = (beta * e_rainfall) + alpha * prev_store
    return max(0.0, store)
end

"""
    calc_flows(prev_quick::Float64, prev_slow::Float64, v_s::Float64,
              e_rainfall::Float64, area::Float64, tau_q::Float64, tau_s::Float64)

Calculate quick and slow flow, and outflow using a linear routing module.

Assumes components are in parallel.

# Arguments
- `prev_quick` : previous quick flow in ML/day
- `prev_slow`  : previous slow flow in ML/day
- `v_s`        : proportional amount that goes to slow flow. v_s <= 1.0
- `e_rainfall` : effective rainfall for t
- `area`      : catchment area
- `tau_q`      : Time constant, quick flow τ value
- `tau_s`      : Time constant slow flow τ value

# Returns
Tuple of (quickflow, slowflow, outflow) in ML/day

# References
- https://wiki.ewater.org.au/display/SD45/IHACRES-CMD+-+SRG
- Croke, B.F.W., Jakeman, A.J. 2004
    A catchment moisture deficit module for the IHACRES rainfall-runoff model,
    Environmental Modelling & Software, 19(1), pp. 1–5.
    doi: 10.1016/j.envsoft.2003.09.001
- Croke, B.F.W., Jakeman, A.J. 2005
    Corrigendum to "A Catchment Moisture Deficit module for the IHACRES
    rainfall-runoff model [Environ. Model. Softw. 19 (1) (2004) 1–5]"
    Environmental Modelling & Software, 20(7), p. 997.
    doi: 10.1016/j.envsoft.2004.11.004
"""
function calc_flows(prev_quick::Float64, prev_slow::Float64, v_s::Float64,
    e_rainfall::Float64, area::Float64, tau_q::Float64, tau_s::Float64
)::Tuple{Float64,Float64,Float64}
    v_q = 1.0 - v_s  # proportional quick flow
    areal_rainfall = e_rainfall * area
    quick = calc_stores(tau_q, prev_quick, v_q, areal_rainfall)
    slow = calc_stores(tau_s, prev_slow, v_s, areal_rainfall)
    outflow = quick + slow
    return (quick, slow, outflow)
end

"""
    routing(gw_vol::Float64, storage_coef::Float64, inflow::Float64, flow::Float64,
            irrig_ext::Float64, gw_exchange::Float64=0.0)

Stream routing, taking into account groundwater interactions and water extractions.

# Arguments
- `gw_vol`       : groundwater store at t-1
- `storage_coef` : groundwater storage factor
- `inflow`       : incoming streamflow (flow from previous node)
- `flow`         : outflow for the node (local flow)
- `irrig_ext`    : volume of irrigation extraction in ML
- `gw_exchange`  : groundwater flux. Defaults to 0.0
                   Negative values represent infiltration into aquifer

# Returns
Tuple of (gw_store, streamflow) in ML/day
"""
function routing(gw_vol::Float64, storage_coef::Float64, inflow::Float64, flow::Float64,
    irrig_ext::Float64, gw_exchange::Float64=0.0)::Tuple{Float64,Float64}
    tmp_gw_store = gw_vol + (inflow + flow + gw_exchange) - irrig_ext

    if tmp_gw_store > 0.0
        # Account for interaction with groundwater system
        c1 = exp(-storage_coef)
        gw_store = c1 * tmp_gw_store
        outflow = (1.0 - c1) * tmp_gw_store
    else
        # Groundwater level is below stream, so no baseflow occurs
        gw_store = tmp_gw_store
        outflow = 0.0
    end

    return (gw_store, outflow)
end

"""
    calc_outflow(flow::Float64, extractions::Float64)

Simple routing modifier to account for extractions.

# Arguments
- `flow`        : sum of quickflow and slowflow in ML/day
- `extractions` : water extractions that occurred in ML/day

# Returns
(flow - extractions), minimum possible is 0.0
"""
function calc_outflow(flow::Float64, extractions::Float64)::Float64
    return max(0.0, flow - extractions)
end

"""
    calc_ft_flows(prev_quick::Float64, prev_slow::Float64, e_rain::Float64,
                 recharge::Float64, area::Float64, a::Float64, b::Float64)

Unit Hydrograph module ported from Fortran.

# Arguments
- `prev_quick` : previous quickflow storage
- `prev_slow`  : previous slowflow storage
- `e_rain`     : effective rainfall in mm
- `recharge`   : recharge amount in mm
- `area`       : catchment area in km²
- `a`          : quickflow storage coefficient, inverse of τ_q such that a := (1/τ_q)
- `b`          : slowflow storage coefficient, inverse of τ_s such that b := (1/τ_s)

# Returns
Tuple of (quick_store, slow_store, outflow)
"""
function calc_ft_flows(prev_quick::Float64, prev_slow::Float64, e_rain::Float64,
    recharge::Float64, area::Float64, a::Float64, b::Float64
)::Tuple{Float64,Float64,Float64}
    tmp_calc = max(0.0, prev_quick + (e_rain * area))

    if tmp_calc > 0.0
        alpha = exp(-a)
        beta = (1.0 - alpha) * tmp_calc
        quick_store = alpha * tmp_calc
        outflow = beta
    else
        quick_store = tmp_calc
        outflow = 0.0
    end

    @assert outflow >= 0.0 "Calculating quick store: Outflow cannot be negative: $outflow"

    slow_store = prev_slow + (recharge * area)
    if slow_store > 0.0
        alpha = exp(-b)
        beta = (1.0 - alpha) * slow_store
        outflow += beta
        slow_store = alpha * slow_store
    end

    @assert outflow >= 0.0 "Calculating slow store: Outflow cannot be negative: $outflow"

    return (quick_store, slow_store, outflow)
end

"""
    calc_ft_level(outflow::Float64, level_params::Vector{Float64})

Calculate stream level from outflow using level parameters.

# Arguments
- `outflow`      : stream outflow
- `level_params` : array of 9 parameters [p1, p2, p3, p4, p5, p6, p7, p8, CTF]

# Returns
Calculated stream level
"""
function calc_ft_level(outflow::Float64, level_params::Vector{Float64})::Float64
    p1, p2, p3, p4, p5, p6, p7, p8, CTF = level_params

    if approx_zero(outflow)
        return CTF  # Return cease to flow level
    end

    level = (exp(p1) * outflow^p2) *
            (1.0 / (
        1.0 + ((outflow / p3)^p4)^(p5 / p4) * exp(p6 / (1.0 + exp(-p7 * p8)) * outflow^p7)
    ))
    level = max(level, 0.0)
    level += CTF  # add Cease to Flow (base height of stream in local datum)

    @assert level >= 0.0 "Stream level cannot be below 0.0, got $level"

    return level
end
