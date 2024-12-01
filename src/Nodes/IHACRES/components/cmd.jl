"""
    calc_cmd(prev_cmd::Float64, rainfall::Float64, et::Float64, effective_rainfall::Float64, recharge::Float64)::Float64

Calculate Catchment Moisture Deficit (CMD).

Min value of CMD is 0.0 and is represented in mm depth.
A value of 0 indicates that the catchment is fully saturated.
A value greater than 0 means that there is a moisture deficit.
"""
function calc_cmd(
    prev_cmd::Float64, rainfall::Float64, et::Float64,
    effective_rainfall::Float64, recharge::Float64
)::Float64
    cmd = prev_cmd + et + effective_rainfall + recharge - rainfall  # units in mm

    return max(0.0, cmd)
end

"""
    calc_linear_interim_cmd(cmd::Float64, param_d::Float64, rainfall::Float64)::Float64

Calculate interim CMD (M_f) in its linear form.

Based on HydroMad implementation and details in references.

# Arguments
- `cmd`: previous Catchment Moisture Deficit (M_{k-1})
- `param_d`: model parameter factor `d`
- `rainfall`: rainfall for current time step in mm

# Returns
interim CMD (M_f)

# References
- Croke, B.F.W., Jakeman, A.J. 2004. A catchment moisture deficit module for the IHACRES rainfall-runoff model,
  Environmental Modelling & Software, 19(1), pp. 1–5. doi: 10.1016/j.envsoft.2003.09.001
- Croke, B.F.W., Jakeman, A.J. 2005. Corrigendum to "A Catchment Moisture Deficit module for the IHACRES
  rainfall-runoff model [Environ. Model. Softw. 19 (1) (2004) 1–5]" Environmental Modelling & Software,
  20(7), p. 997. doi: 10.1016/j.envsoft.2004.11.004
"""
function calc_linear_interim_cmd(cmd::Float64, param_d::Float64, rainfall::Float64)::Float64
    if cmd < param_d
        return cmd * exp(-rainfall / param_d)
    elseif cmd < (param_d + rainfall)
        return param_d * exp((-rainfall + cmd - param_d) / param_d)
    end

    return cmd - rainfall
end

"""
    calc_trig_interim_cmd(cmd::Float64, param_d::Float64, rainfall::Float64)::Float64

Calculate interim CMD (M_f) in its trigonometric form.

Based on HydroMad implementation and details in references.

# Arguments
- `cmd`: previous Catchment Moisture Deficit (M_{k-1})
- `param_d`: model parameter `d`
- `rainfall`: rainfall for current time step in mm

# Returns
- interim CMD (M_f)
"""
function calc_trig_interim_cmd(cmd::Float64, param_d::Float64, rainfall::Float64)::Float64
    if cmd < param_d
        Mf = 1.0 / tan((cmd / param_d) * (π / 2.0))
        return (2.0 * param_d / π) * atan(1.0 / (π * rainfall / (2.0 * param_d) + Mf))
    elseif rainfall < (param_d + rainfall)
        return (2.0 * param_d / π) * atan(2.0 * param_d / (π * (param_d - cmd + rainfall)))
    else
        return cmd - rainfall
    end
end

"""
    calc_ft_interim_cmd(cmd::Float64, rain::Float64, d::Float64, d2::Float64, alpha::Float64)

Direct port of original Fortran implementation to calculate interim CMD (M_f).
Calculates estimates of effective rainfall and recharge as a by-product.

# Arguments
- `cmd`: previous Catchment Moisture Deficit (M_{k-1})
- `rain`: rainfall for time step in mm
- `d`: flow threshold value
- `d2`: scaling factor applied to `d`
- `alpha`: scaling factor

# Returns
A tuple of (interim CMD value, effective rainfall, recharge)
"""
function calc_ft_interim_cmd(cmd::Float64, rain::Float64, d::Float64, d2_factor::Float64, alpha::Float64)::Tuple{Float64, Float64, Float64}
    # Initialize variables
    d2 = d * d2_factor
    tmp_cmd = cmd
    e_rain = 0.0
    recharge = 0.0

    # Early return if no rain
    if isapprox(rain, 0.0)
        return (cmd, e_rain, recharge)
    end

    # Calculate Mf (new CMD value)
    if tmp_cmd > (d2 + rain)
        # CMD never reaches d2, so all rain is effective
        Mf = tmp_cmd - rain
    else
        # Handle more complex cases
        if tmp_cmd > d2
            tmp_rain = rain - (tmp_cmd - d2)  # leftover rain after reaching d2 threshold
            tmp_cmd = d2
        else
            tmp_rain = rain
        end

        # Calculate threshold value d1a
        d1a = d * (2.0 - exp(-(rain / 50.0)^2))

        if tmp_cmd > d1a
            eps = d2 / (1.0 - alpha)

            # Amount of rain necessary to get to threshold d
            depth_to_d = eps * log((alpha + tmp_cmd / eps) / (alpha + d1a / eps))

            if depth_to_d >= tmp_rain
                lam = exp(tmp_rain * (1.0 - alpha) / d2)  # lambda
                epsilon = alpha * eps

                Mf = tmp_cmd / lam - epsilon * (1.0 - 1.0 / lam)
                e_rain = 0.0
            else
                if tmp_cmd > d1a
                    tmp_rain = tmp_rain - depth_to_d
                end

                tmp_cmd = d1a
                gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
                Mf = tmp_cmd * exp(-tmp_rain * gamma)
                e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (Mf - tmp_cmd))
            end
        else
            gamma = (alpha * d2 + (1.0 - alpha) * d1a) / (d1a * d2)
            Mf = tmp_cmd * exp(-tmp_rain * gamma)
            e_rain = alpha * (tmp_rain + 1.0 / d1a / gamma * (Mf - tmp_cmd))
        end

        recharge = rain - (cmd - Mf) - e_rain
    end

    return (Mf, e_rain, recharge)
end
