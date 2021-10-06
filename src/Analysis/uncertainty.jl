using Dates
using Statistics, DataFrames, Bootstrap
import ..ME


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

        Q1, Q2, Q3 = quantile(obs_gi, [0.25, 0.5, 0.75])
        IQR_lim = 1.5*(Q3 - Q1)
        min_R = Q1 - IQR_lim
        max_R = Q3 + IQR_lim

        min_section[i] = min_R
        max_section[i] = max_R
        x_section[i] = Q2
    end

    whisker_range = round(mean(max_section .- min_section), digits=2)

    return x_section, min_section, max_section, whisker_range
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

    x_section, min_section, max_section, whisker_range = temporal_uncertainty(dates, obs, period)

    # Increase weight on observations with lower uncertainty
    weights = 1.0 .- (abs.(x_section) ./ (threshold * std(x_section))).^2
    weights = max.(weights, 0.0)

    # Increase weight on observations that are more uncertain
    # abs_xsect = abs.(x_section)
    # weights = abs_xsect ./ (abs_xsect .+ (threshold * std(x_section)))

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

    return x_section, min_section, max_section, whisker_range, weighted
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
"""
function temporal_uncertainty(dates, obs, sim, period::Function, func::Function)
    df = DataFrame(Date=dates, Observed=obs, Modeled=sim)
    sp = sort(unique(period.(dates)))

    # Remove leap days (when using monthday)
    deleteat!(sp, findall(x -> x == (2,29), sp))

    x_section = fill(0.0, length(sp))
    min_section = fill(0.0, length(sp))
    max_section = fill(0.0, length(sp))
    for (i, obs_i) in enumerate(sp)
        obs_g = df[in([obs_i]).(period.(df.Date)), :]
        obs_gi = obs_g.Observed
        sim_gi = obs_g.Modeled

        tmp = func.([[x] for x in obs_gi], [[x] for x in sim_gi])

        Q1, Q2, Q3 = quantile(tmp, [0.25, 0.5, 0.75])
        IQR_lim = 1.5*(Q3 - Q1)
        min_R = Q1 - IQR_lim
        max_R = Q3 + IQR_lim

        min_section[i] = min_R
        max_section[i] = max_R
        x_section[i] = Q2
    end

    whisker_range = round(mean(max_section .- min_section), digits=2)

    return x_section, min_section, max_section, whisker_range
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
function temporal_uncertainty(dates, obs, sim; period::Function=monthday, func::Function=ME, threshold::Float64=1.5)

    x_section, min_section, max_section, whisker_range = temporal_uncertainty(dates, obs, sim, period, func)

    # Increase weight on observations with lower uncertainty
    weights = 1.0 .- (abs.(x_section) ./ (threshold * std(x_section))).^2
    weights = max.(weights, 0.0)

    # Increase weight on observations that are more uncertain
    # abs_xsect = abs.(x_section)
    # weights = abs_xsect ./ (abs_xsect .+ (threshold * std(x_section)))

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

    return x_section, min_section, max_section, whisker_range, weighted
end