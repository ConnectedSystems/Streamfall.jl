using Dates
using Statistics, DataFrames, Bootstrap
import ..ME, ..MAE


"""
Conduct temporal uncertainty analysis on observations.

Calculates the standard boxplot values for given subperiods (median, -1.5*IQR, +1.5*IQR).

Notes: Ignores leap days when using monthday subperiods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : Q1 - 1.5*IQR of metric value for period
     - max_section : Q3 + 1.5*IQR of metric value for period
     - whisker_range : Mean range indicated by (max_section - min_section)
"""
function temporal_uncertainty(dates, obs, period::Function)
    df = DataFrame(Date=dates, Observed=obs)
    sp = sort(unique(period.(dates)))

    # Remove leap days (when using monthday)
    deleteat!(sp, findall(x -> x == (2,29), sp))

    x_section = fill(0.0, length(sp))
    min_section = fill(0.0, length(sp))
    max_section = fill(0.0, length(sp))
    mad_section = fill(0.0, length(sp))
    for (i, obs_i) in enumerate(sp)
        obs_g = df[in([obs_i]).(period.(df.Date)), :]
        obs_gi = obs_g.Observed

        ab_min, lower_q, upper_q, ab_max = quantile(obs_gi, [0.0, 0.05, 0.95, 1.0])
        low_section[i] = lower_q
        upp_section[i] = upper_q
        min_section[i] = ab_min
        max_section[i] = ab_max
        mae_section[i] = MAE(obs_gi, sim_gi)
        x_section[i] = mean(obs_gi)

        mad_section[i] = mean(abs.(obs_gi .- mean(obs_gi)))
    end

    whisker_range = round(mean(max_section .- min_section), digits=2)

    # diff_x = diff(x_section)
    # normed_diff = (diff_x .- mean(diff_x)) ./ std(diff_x)
    # roughness = mean(diff(normed_diff).^2 ./ 4)  # 0 = smooth, 1 = maximal roughness
    # mean_mad = mean(mad_section)

    cv = std(whisker_range) / mean(whisker_range)
    return x_section, low_section, upp_section, whisker_range, cv, roughness
end


"""
Conduct temporal uncertainty analysis on observations.

Calculates the standard boxplot values for given subperiods (median, -1.5*IQR, +1.5*IQR).

Periods of high uncertainty, defined as 1.5 std from mean tolerance interval by default,
are given weights > 1.0, increasing with distance. Those within this interval is given lesser weight.

Notes: Ignores leap days when using monthday periods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : Q1 - 1.5*IQR of metric value for period
     - max_section : Q3 + 1.5*IQR of metric value for period
     - whisker_range : Mean range indicated by max_section - min_section
     - weighted : weights (0 to 1) indicating relative to `threshold * std(x_section)`
"""
function temporal_uncertainty(dates, obs; period::Function=monthday, threshold::Float64=1.5)

    x_section, low_section, upp_section, whisker_range, roughness = temporal_uncertainty(dates, obs, period)

    # Increase weight on observations with lower uncertainty
    # weights = 1.0 .- (abs.(x_section) ./ (threshold * std(x_section))).^2
    # weights = max.(weights, 0.0)

    # Increase weight on observations that are more uncertain
    abs_xsect = abs.(x_section)
    weights = abs_xsect ./ (abs_xsect .+ (threshold * std(x_section)))

    weighted = fill(1.0, length(obs))
    doy = Dates.dayofyear.(dates)
    mds = Dates.monthday.(dates)
    for (idx, (d, md)) in enumerate(zip(doy, mds))
        if md == (2, 29)
            weighted[idx] = weights[d+1]
            continue
        elseif d == 366
            weighted[idx] = weights[d-1]
            continue
        end

        weighted[idx] = weights[d]
    end

    return x_section, low_section, upp_section, whisker_range, weighted, roughness
end


"""
Conduct temporal uncertainty analysis on simulated results relative to observations
using the provided function `func`.

Calculates the standard boxplot values for given subperiods (median, -1.5*IQR, +1.5*IQR).

Notes: Ignores leap days when using monthday subperiods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : Q1 - 1.5*IQR of metric value for period
     - max_section : Q3 + 1.5*IQR of metric value for period
     - whisker_range : Mean range indicated by (max_section - min_section)
     - roughness : Value 0 to 1 where 0 is smooth and 1 is maximal roughness (lower values better)
"""
function temporal_uncertainty(dates, obs, sim, period::Function, func::Function)
    df = DataFrame(Date=dates, Observed=obs, Modeled=sim)
    sp = sort(unique(period.(dates)))

    # Remove leap days (when using monthday)
    deleteat!(sp, findall(x -> x == (2,29), sp))

    x_section = fill(0.0, length(sp))
    min_section = fill(0.0, length(sp))
    max_section = fill(0.0, length(sp))
    low_section = fill(0.0, length(sp))
    upp_section = fill(0.0, length(sp))
    mae_section = fill(0.0, length(sp))
    for (i, obs_i) in enumerate(sp)
        obs_g = df[in([obs_i]).(period.(df.Date)), :]
        obs_gi = obs_g.Observed
        sim_gi = obs_g.Modeled

        tmp = func.([[x] for x in obs_gi], [[x] for x in sim_gi])

        ab_min, lower_q, upper_q, ab_max = quantile(tmp, [0.0, 0.05, 0.95, 1.0])

        low_section[i] = lower_q
        upp_section[i] = upper_q
        min_section[i] = ab_min
        max_section[i] = ab_max
        mae_section[i] = MAE(obs_gi, sim_gi)
        x_section[i] = mean(tmp)
    end

    whisker_range = round(mean(max_section .- min_section), digits=2)

    # diff_x = diff(x_section)
    # normed_diff = (diff_x .- mean(diff_x)) ./ std(diff_x)
    # roughness = mean(diff(normed_diff).^2 ./ 4)  # 0 = smooth, 1 = maximal roughness
    cv = std(whisker_range) / mean(whisker_range)

    return x_section, low_section, upp_section, whisker_range, cv
end


"""
Conduct temporal uncertainty analysis on observations.

Calculates the standard boxplot values for given subperiods (median, -1.5*IQR, +1.5*IQR).

Periods of high uncertainty, defined as 1.5 std from mean tolerance interval by default,
are given weights > 1.0, increasing with distance. Those within this interval is given lesser weight.

Notes: Ignores leap days when using monthday periods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : Q1 - 1.5*IQR of metric value for period
     - max_section : Q3 + 1.5*IQR of metric value for period
     - whisker_range : Mean range indicated by max_section - min_section
     - weighted : weights (0 to 1) relative to `threshold * std(x_section)`
"""
function temporal_uncertainty(dates, obs, sim; period::Function=monthday, func::Function=ME, threshold::Float64=1.5)

    x_section, min_section, max_section, whisker_range, roughness = temporal_uncertainty(dates, obs, sim, period, func)

    # Increase weight on observations with lower uncertainty
    # weights = 1.0 .- (abs.(x_section) ./ (threshold * std(x_section))).^2
    # weights = max.(weights, 0.0)

    # Increase weight on observations with larger error bounds
    # abs_xsect = abs.(x_section)
    # weights = abs_xsect ./ (abs_xsect .+ (threshold * std(x_section)))

    # inverted standardized
    weights = 1.0 / (x_section / (threshold * std(x_section)))

    weighted = fill(1.0, length(obs))
    doy = Dates.dayofyear.(dates)
    mds = Dates.monthday.(dates)
    for (idx, (d, md)) in enumerate(zip(doy, mds))
        if md == (2, 29)
            weighted[idx] = weights[d+1]
            continue
        elseif d == 366
            weighted[idx] = weights[d-1]
            continue
        end

        weighted[idx] = weights[d]
    end

    return x_section, min_section, max_section, whisker_range, roughness, weighted
end