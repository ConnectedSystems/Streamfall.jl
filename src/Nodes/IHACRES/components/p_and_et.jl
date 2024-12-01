export calc_effective_rainfall, calc_ET_from_E, calc_ET, calc_ET_from_T

"""
    calc_effective_rainfall(rainfall::Float64, cmd::Float64, d::Float64, d2::Float64, n::Float64=0.1)

Estimate effective rainfall.

# References
- Croke, B.F.W., Jakeman, A.J. 2004
  A catchment moisture deficit module for the IHACRES rainfall-runoff model,
  Environmental Modelling & Software, 19(1), pp. 1–5.
  doi: 10.1016/j.envsoft.2003.09.001

- Croke, B.F.W., Jakeman, A.J. 2005
  Corrigendum to "A Catchment Moisture Deficit module for the IHACRES
  rainfall-runoff model [Environ. Model. Softw. 19 (1) (2004) 1–5]"
  Environmental Modelling & Software, 20(7), p. 997.
  doi: https://doi.org/10.1016/j.envsoft.2004.11.004

# Arguments
- `rainfall`: rainfall for time step
- `cmd`: previous CMD value
- `d`: threshold value
- `d2`: scaling factor applied to `d`
- `n`: scaling factor (default = 0.1). Default value is suitable for most cases (Croke & Jakeman, 2004)

# Returns
- Effective rainfall value
"""
function calc_effective_rainfall(rainfall::Float64, cmd::Float64, d::Float64,
                               d2::Float64, n::Float64=0.1)::Float64
    scaled_d = d * d2

    if cmd > scaled_d
        e_rainfall = rainfall
    else
        f1 = min(1.0, cmd / d)
        f2 = min(1.0, cmd / scaled_d)
        e_rainfall = rainfall * ((1.0 - n) * (1.0 - f1) + (n * (1.0 - f2)))
    end

    return max(0.0, e_rainfall)
end

"""
    calc_ET_from_E(e::Float64, evap::Float64, Mf::Float64, f::Float64, d::Float64)

Calculate evapotranspiration from evaporation.

# Arguments
- `e`: temperature to PET conversion factor (a stress threshold)
- `evap`: evaporation for given time step
- `Mf`: Catchment Moisture Deficit prior to accounting for ET losses (`Mf`)
- `f`: calibrated parameter that acts as a multiplication factor on `d`
- `d`: flow threshold factor

# Returns
- Estimate of ET
"""
function calc_ET_from_E(e::Float64, evap::Float64, Mf::Float64, f::Float64,
                       d::Float64)::Float64
    param_g = f * d
    et = e * evap

    if Mf > param_g
        et *= min(1.0, exp((1.0 - Mf/param_g)*2.0))
    end

    return max(0.0, et)
end

# """
#     calc_ET(e::Float64, evap::Float64, Mf::Float64, f::Float64, d::Float64)

# !!! warning
#     Deprecated function - use [`calc_ET_from_E`](@ref) instead.
# """
# function calc_ET(e::Float64, evap::Float64, Mf::Float64, f::Float64, d::Float64)::Float64
#     return calc_ET_from_E(e, evap, Mf, f, d)
# end

"""
    calc_ET_from_T(e::Float64, T::Float64, Mf::Float64, f::Float64, d::Float64)

Calculate evapotranspiration based on temperature data.

Parameters `f` and `d` are used to calculate `g`, the value of the CMD
which the ET rate will begin to decline due to insufficient
water availability for plant transpiration.

# Arguments
- `e`: temperature to PET conversion factor (a stress threshold)
- `T`: temperature in degrees C
- `Mf`: Catchment Moisture Deficit prior to accounting for ET losses (`Mf`)
- `f`: multiplication factor on `d`
- `d`: flow threshold factor

# Returns
- Estimate of ET from temperature (for catchment area)
"""
function calc_ET_from_T(e::Float64, T::Float64, Mf::Float64, f::Float64,
                       d::Float64)::Float64
    # temperature can be negative, so we have a min cap of 0.0
    if T <= 0.0
        return 0.0
    end

    param_g = f * d
    et = e * T * min(1.0, exp(2.0 * (1.0 - (Mf / param_g))))

    return max(0.0, et)
end
