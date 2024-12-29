using Plots, StatsPlots
using Plots.Measures
using DataFrames, Dates, Statistics, Distributions, LaTeXStrings
import Bootstrap: bootstrap, BalancedSampling

import .Analysis: TemporalCrossSection


function quickplot(node::NetworkNode)
    fig = plot(node.outflow)
    return fig
end


function quickplot(node::NetworkNode, climate::Climate)
    date = timesteps(climate)

    @assert length(date) == length(node.outflow) || "Date length and result lengths do not match!"

    fig = plot(date, node.outflow)

    return fig
end


function quickplot(obs, node::NetworkNode, climate::Climate, label="", log=true; burn_in=1, limit=nothing, metric=Streamfall.mKGE)
    return quickplot(obs, node.outflow, climate, label, log; burn_in=burn_in, limit=limit, metric=metric)
end

function quickplot(obs::DataFrame, sim::Vector, climate::Climate, label="", log=true; burn_in=1, limit=nothing, metric=Streamfall.mKGE)
    return quickplot(Matrix(obs[:, Not("Date")])[:, 1], sim, climate, label, log; burn_in, limit, metric)
end
function quickplot(obs::Vector, sim::Vector, climate::Climate, label="", log=true; burn_in=1, limit=nothing, metric=Streamfall.mKGE)
    date = timesteps(climate)
    last_e = !isnothing(limit) ? limit : lastindex(obs)
    show_range = burn_in:last_e
    return quickplot(obs[show_range], sim[show_range], date[show_range], label, log; metric=metric)
end


function quickplot(obs::Vector, sim::Vector, xticklabels::Vector, label="Modeled", log=true; metric=Streamfall.mKGE)
    @assert length(xticklabels) == length(obs) || "x-axis tick label length and observed lengths do not match!"
    @assert length(xticklabels) == length(sim) || "x-axis tick label length and simulated lengths do not match!"

    score = round(metric(obs, sim), digits=4)
    metric_name = String(Symbol(metric))

    if log
        # Add small constant in case of 0-flow
        obs = obs .+ 1e-2
        sim = sim .+ 1e-2
    end

    label = "$(label) ($(metric_name): $(score))"
    fig = plot(xticklabels, obs,
                label="Observed",
                legend=:best,
                ylabel="Streamflow",
                xlabel="Date",
                fg_legend=:transparent,
                bg_legend=:transparent)
    plot!(xticklabels, sim, label=label, alpha=0.7)

    if log
        # modify yaxis
        yaxis!(fig, :log10)
    end

    qqfig = qqplot(obs, sim, legend=false, markerstrokewidth=0, alpha=0.7, xlabel="Observed", ylabel="Modeled")

    if log
        xaxis!(qqfig, :log10)
        yaxis!(qqfig, :log10)
    end

    combined = plot(fig, qqfig, size=(1000, 500), left_margin=5mm, bottom_margin=5mm, layout=(1,2))

    return combined
end


"""
    plot_residuals(obs::Array, sim::Array; xlabel="", ylabel="", title="")

Plot residual between two sequences.

# Arguments
- x : x-axis data
- y : y-axis data
- xlabel : x-axis label
- ylabel : y-axis label
- title : title text
"""
function plot_residuals(x::Array, y::Array; xlabel="", ylabel="", title="")
    # 1:1 Plot
    fig_1to1 = scatter(x, y, legend=false,
                       markerstrokewidth=0, markerstrokealpha=0, alpha=0.2)
    plot!(x, y, color=:red, markersize=0.1, markerstrokewidth=0,
    xlabel=xlabel, ylabel=ylabel, title=title)

    return fig_1to1
end


"""Symmetrical log values.

https://kar.kent.ac.uk/32810/2/2012_Bi-symmetric-log-transformation_v5.pdf
https://discourse.julialang.org/t/symmetrical-log-plot/45709/3
"""
function symlog(y)
    return sign.(y) .* log10.(1.0 .+ abs.(y))
end


"""
    temporal_cross_section(dates, obs; ylabel=nothing, period::Function=monthday)

Provides indication of temporal variation and uncertainty across time, grouped by `period`.

Notes:
Assumes daily data.
Filters out leap days.

# Arguments
- dates : Date of each observation
- obs : observed data
- ylabel : Optional replacement ylabel. Uses name of `func` if not provided.
- `period::Function` : Method from `Dates` package to group (defaults to `monthday`)
"""
function temporal_cross_section(dates, obs;
                                title="", ylabel="ME", label=nothing,
                                period::Function=monthday,
                                kwargs...)  # show_extremes::Bool=false,
    if isnothing(label)
        label = ylabel
    end

    arg_keys = keys(kwargs)
    format_func = y -> y
    logscale = [:log, :log10]
    tmp = nothing

    xsect_res = TemporalCrossSection(dates, obs, period)

    if :yscale in arg_keys || :yaxis in arg_keys
        tmp = (:yscale in arg_keys) ? kwargs[:yscale] : kwargs[:yaxis]

        if tmp in logscale
            log_obs = symlog(copy(obs))

            # Format function for y-axis tick labels (e.g., 10^x)
            format_func = y -> (y != 0) ? L"%$(Int(round(sign(y)) * 10))^{%$(round(abs(y), digits=1))}" : L"0"

            log_xsect_res = TemporalCrossSection(dates, log_obs, period)
            target = log_xsect_res.cross_section
        else
            target = xsect_res.cross_section
        end
    else
        target = xsect_res.cross_section
    end

    x_section = target.median
    lower_75 = target.lower_75
    upper_75 = target.upper_75
    lower_95 = target.lower_95
    upper_95 = target.upper_95

    sp = target.subperiod
    xlabels = join.(sp, "-")

    if !isnothing(tmp) & (tmp in logscale)
        # Remove keys so later plotting methods do not error
        kwargs = Dict(kwargs)
        delete!(kwargs, :yscale)
        delete!(kwargs, :yaxis)
    end

    # Display indicator values using original data instead of log-transformed data
    med = xsect_res.cross_section.median
    m_ind = round(mean(med), digits=2)
    sd_ind = round(std(med), digits=2)

    wr75_m_ind = round(xsect_res.mean_75, digits=2)
    wr75_sd_ind = round(xsect_res.std_75, digits=2)
    wr95_m_ind = round(xsect_res.mean_95, digits=2)
    wr95_sd_ind = round(xsect_res.std_95, digits=2)

    fig = plot(xlabels, lower_95, fillrange=upper_95, color="lightblue", alpha=0.3, label="CI₉₅ μ: $(wr95_m_ind), σ: $(wr95_sd_ind)", linealpha=0)
    plot!(fig, xlabels, lower_75, fillrange=upper_75, color="lightblue", alpha=0.5, label="CI₇₅ μ: $(wr75_m_ind), σ: $(wr75_sd_ind)", linealpha=0)
    plot!(fig, xlabels, x_section,
            label="Mean of $(label) μ: $(m_ind), σ: $(sd_ind)",
            color="black",
            xlabel=nameof(period),
            ylabel=ylabel,
            legend=:bottomleft,
            legendfont=Plots.font(10),
            fg_legend=:transparent,
            bg_legend=:transparent,
            left_margin=5mm,
            bottom_margin=5mm,
            title=title,
            yformatter=format_func;
            kwargs...)

    # if show_extremes
    #     scatter!(fig, xlabels, min_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
    #     scatter!(fig, xlabels, max_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
    # end

    return fig
end


"""
    temporal_cross_section(
        dates, obs, sim;
        title="", ylabel="Median Error", label=nothing, period::Function=monthday, kwargs...
    )

Provides indication of temporal variation and uncertainty across time, grouped by `period`.

Notes:
Assumes daily data.
Filters out leap days.

# Arguments
- `dates` : Date of each observation
- `obs` : Observed data
- `sim` : Simulated data
- `title` : Optional plot title. Blank if not provided.
- `ylabel` : Optional replacement ylabel. Uses name of `period` if not provided.
- `label` : Optional legend label. Uses `ylabel` if not provided.
- `period` : Method from `Dates` package to group (defaults to `month`)
"""
function temporal_cross_section(
    dates, obs, sim;
    title="", ylabel="Median Error", label=nothing, period::Function=monthday, kwargs...
)  # show_extremes::Bool=false,
    if isnothing(label)
        label = ylabel
    end

    arg_keys = keys(kwargs)
    format_func = y -> y
    logscale = [:log, :log10]
    tmp = nothing

    xsect_res = TemporalCrossSection(dates, obs, sim, period)
    target = xsect_res.cross_section

    if :yscale in arg_keys || :yaxis in arg_keys
        tmp = (:yscale in arg_keys) ? kwargs[:yscale] : kwargs[:yaxis]

        if tmp in logscale
            log_obs = symlog(copy(obs))
            log_sim = symlog(copy(sim))

            # Format function for y-axis tick labels (e.g., 10^x)
            format_func = y -> (y != 0) ? L"%$(Int(round(sign(y)) * 10))^{%$(round(abs(y), digits=1))}" : L"0"

            log_xsect_res = TemporalCrossSection(dates, log_obs, log_sim, period)
            target = log_xsect_res.cross_section
        end
    end

    x_section = target.median
    lower_75 = target.lower_75
    upper_75 = target.upper_75
    lower_95 = target.lower_95
    upper_95 = target.upper_95

    sp = target.subperiod
    xlabels = join.(sp, "-")

    if !isnothing(tmp) & (tmp in logscale)
        # Remove keys so later plotting methods do not error
        kwargs = Dict(kwargs)
        delete!(kwargs, :yscale)
        delete!(kwargs, :yaxis)
    end

    xsect_df = xsect_res.cross_section

    # Display indicator values using original data instead of log-transformed data
    med = xsect_df.median
    m_ind = round(mean(med), digits=2)
    sd_ind = round(std(med), digits=2)

    wr75_m_ind = round(xsect_res.mean_75, digits=2)
    wr75_sd_ind = round(xsect_res.std_75, digits=2)
    wr95_m_ind = round(xsect_res.mean_95, digits=2)
    wr95_sd_ind = round(xsect_res.std_95, digits=2)

    fig = plot(xlabels, lower_95, fillrange=upper_95, color="lightblue", alpha=0.3, label="CI₉₅ μ: $(wr95_m_ind), σ: $(wr95_sd_ind)", linealpha=0)
    plot!(fig, xlabels, lower_75, fillrange=upper_75, color="lightblue", alpha=0.5, label="CI₇₅ μ: $(wr75_m_ind), σ: $(wr75_sd_ind)", linealpha=0)
    plot!(
        fig, xlabels, x_section,
        label="$(label) μ: $(m_ind), σ: $(sd_ind)",
        color="black",
        xlabel=nameof(period),
        ylabel=ylabel,
        legend=:bottomleft,
        legendfont=Plots.font(10),
        fg_legend=:transparent,
        bg_legend=:transparent,
        left_margin=5mm,
        bottom_margin=5mm,
        title=title,
        yformatter=format_func;
        kwargs...
    )

    # if show_extremes
    #     scatter!(fig, xlabels, min_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
    #     scatter!(fig, xlabels, max_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
    # end

    return fig
end
