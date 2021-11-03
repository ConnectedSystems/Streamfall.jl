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
    for (i, obs_i) in enumerate(sp)
        obs_g = df[in([obs_i]).(period.(df.Date)), :]
        obs_gi = obs_g.Observed

        ab_min, lower_q, upper_q, ab_max = quantile(obs_gi, [0.0, 0.05, 0.95, 1.0])
        low_section[i] = lower_q
        upp_section[i] = upper_q
        min_section[i] = ab_min
        max_section[i] = ab_max
        # mae_section[i] = MAE(obs_gi, sim_gi)
        x_section[i] = median(obs_gi)
    end

    whisker_range = round(mean(upp_section .- low_section), digits=2)

    # diff_x = diff(x_section)
    # normed_diff = (diff_x .- mean(diff_x)) ./ std(diff_x)
    # roughness = mean(diff(normed_diff).^2 ./ 4)  # 0 = smooth, 1 = maximal roughness
    # mean_mad = mean(mad_section)

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
     - weighted : weights (0 to 1) indicating relative to `threshold * std(x_section)`
"""
function temporal_uncertainty(dates, obs; period::Function=monthday, min_weight::Float64=0.3)

    x_section, low_section, upp_section, whisker_range, cv_r = temporal_uncertainty(dates, obs, period)

    # use min-max scaling to indicate where to place greater weight
    weights = 1 .+ (whisker_range .- minimum(whisker_range)) ./ (maximum(whisker_range) - minimum(whisker_range))
    weights[weights .< (1.0 + min_weight)] .= (1.0 + min_weight)

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

    return x_section, low_section, upp_section, whisker_range, weighted, cv_r
end


"""
Conduct temporal uncertainty analysis on simulated results relative to observations
using the provided function `func`.

Calculates the standard boxplot values for given subperiods (median, -1.5*IQR, +1.5*IQR).

Notes: Ignores leap days when using monthday subperiods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : minimum score for given period
     - max_section : maximum score for given period
     - whisker_range : Mean range indicated by q0.95 - q0.05 (i.e., 95th - 5th percentile)
     - cv_r : coefficient of variation for whisker range. Indicates variability of 
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
    # mae_section = fill(0.0, length(sp))
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
        # mae_section[i] = MAE(obs_gi, sim_gi)
        x_section[i] = median(tmp)
    end

    whisker_range = upp_section .- low_section

    # diff_x = diff(x_section)
    # normed_diff = (diff_x .- mean(diff_x)) ./ std(diff_x)
    # roughness = mean(diff(normed_diff).^2 ./ 4)  # 0 = smooth, 1 = maximal roughness
    cv_r = std(whisker_range) / mean(whisker_range)

    return x_section, low_section, upp_section, min_section, max_section, whisker_range, cv_r, stderror(x_section)
end


"""
Conduct temporal uncertainty analysis on observations.

Calculates boxplot-like values for given subperiods (mean, q0.05, q0.95).

A temporal cross-section is taken for each temporal period (default: month-day)
and greater weights are placed on periods of greater uncertainty.

Min-max scaled range of (q0.95 - q0.05) indicates where greater weight (``w``) is placed,
such that ``w=1+α`` for the most confident time period (i.e., lowest range) and ``w=2``
for the least confident time period, where ``α`` is the minimum additional weight to apply
(defaults to 0).

Notes: Ignores leap days when using monthday periods.

# Returns
tuple: 
     - x_section : Median of observations for a given period (default: day of year)
     - min_section : minimum score for a given period
     - max_section : maximum score for a given period
     - whisker_range : Range indicated by q0.95 - q0.05 (95th - 5th quantile)
     - cv_r : coefficient of variation across whisker_range
     - weighted : weights (1 to 2), with minimum weight of 1 + `min_weight`
"""
function temporal_uncertainty(dates, obs, sim; period::Function=monthday, func::Function=ME, min_weight::Float64=1.0)

    x_section, low_section, upp_section, min_section, max_section, whisker_range, cv_r, std_error = temporal_uncertainty(dates, obs, sim, period, func)

    # use min-max scaling to indicate where to place greater weight
    max_range = max_section .- min_section
    min_of_range = minimum(max_range)
    weights = 1 .+ ((max_range .- min_of_range) ./ (maximum(max_range) - min_of_range))
    weights[weights .< (1.0 + min_weight)] .= (1.0 + min_weight)

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

    return x_section, low_section, upp_section, min_section, max_section, whisker_range, cv_r, std_error, weighted
end