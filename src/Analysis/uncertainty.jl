using Dates
using Statistics, DataFrames, Bootstrap
using OrderedCollections
import ..ME, ..MAE


mutable struct TemporalCrossSection{A<:Dates.Date, B<:Real, D<:OrderedDict}
    dates::Array{A}
    ts::Array{B}
    subperiods::D
    cross_section::Union{DataFrame, Nothing}
end


function TemporalCrossSection(dates::Array, ts::Array, period::Function=monthday)
    df = DataFrame(Date=dates, data=ts)
    sp = sort(unique(period.(dates)))

    # Remove leap days (when using monthday)
    deleteat!(sp, findall(x -> x == (2,29), sp))

    res = OrderedDict(
        sp_i => df[period.(df.Date) .== [sp_i], :].data for sp_i in sp
    )

    return TemporalCrossSection(dates, ts, res, nothing)
end


"""
TemporalCrossSection constructor that handles two time series
"""
function TemporalCrossSection(dates::Array, obs::Array, sim::Array,
                              period::Function=monthday)
    Y = ME.(obs, sim)
    return TemporalCrossSection(dates, Y, period)
end


function cross_section(tr::TemporalCrossSection)
    sp = tr.subperiods

    # Array holding boxplot-like data
    # No. of Dims: (75, 95) * (2) + median + min/max = 7
    n_cols = 9  # min, -q95, -q75, q50, q75, q95, max, mean, std
    boxframe = fill(0.0, length(keys(sp)), n_cols)

    target_pts = [0.0, 0.025, 0.125, 0.875, 0.975, 1.0]
    for (i, (ki, vi)) in enumerate(sp)
        ab_min, lower_95, lower_75, upper_75, upper_95, ab_max = quantile(vi, target_pts)

        boxframe[i, 1:3] .= ab_min, lower_95, lower_75
        boxframe[i, 4] = median(vi)
        boxframe[i, 5:7] .= upper_75, upper_95, ab_max
        boxframe[i, 8:9] .= mean(vi), std(vi)
    end

    cols = [:abs_min, :lower_95, :lower_75, :median, 
            :upper_75, :upper_95, :abs_max, :mean, :std]

    x_section = DataFrame(boxframe, cols)
    x_section[:, :subperiod] = collect(keys(sp))

    tr.cross_section = x_section

    return x_section
end


"""
    offsets(tr::TemporalCrossSection, period::Function=monthday)

Determine (median) offsets, assuming mean error-based temporal cross-sectional analysis 
was performed. Only supports daily time series.
"""
function offsets(tr::TemporalCrossSection)
    x_section = tr.cross_section
    offset = x_section.median
    dates = tr.dates

    offset_vals = zeros(length(dates))

    doy = Dates.dayofyear.(dates)
    mds = monthday.(dates)
    for (idx, (d, md)) in enumerate(zip(doy, mds))
        if md == (2, 29) || d == 366
            offset_vals[idx] = offset[d-1] * -1.0
            continue
        end

        offset_vals[idx] = offset[d] * -1.0
    end

    return offset_vals
end


function Base.getproperty(a::TemporalCrossSection, v::Symbol)
    if v == :cross_section
        if isnothing(getfield(a, :cross_section))
            boxframe = cross_section(a)
            return boxframe
        end
    end

    if v == :offsets
        return offsets(a)
    end

    if in(v, [:cv_median, :cv_95, :cv_75, :mean_95, :mean_75, :std_75, :std_95])
        boxframe = a.cross_section

        if v == :cv_median
            med = boxframe[:, :median]
            cv_med = std(med) ./ mean(med)
            return cv_med
        elseif in(v, [:cv_95, :mean_95, :std_95])
            upper = boxframe[:, :upper_95]
            lower = boxframe[:, :lower_95]
            wr = upper .- lower

            if v == :cv_95
                return std(wr) ./ mean(wr)
            elseif v == :mean_95
                return mean(wr)
            elseif v == :std_95
                return std(wr, corrected=false)
            end
        elseif in(v, [:cv_75, :mean_75, :std_75])
            upper = boxframe[:, :upper_75]
            lower = boxframe[:, :lower_75]
            wr = upper .- lower

            if v == :cv_75
                return std(wr) ./ mean(wr)
            elseif v == :mean_75
                return mean(wr)
            elseif v == :std_75
                return std(wr, corrected=false)
            end
        end
    else
        return getfield(a, v)
    end
end
