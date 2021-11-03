using Plots, StatsPlots
using Plots.Measures
using DataFrames, Dates, Statistics, Distributions, LaTeXStrings
import StatsBase: ecdf
import Bootstrap: bootstrap, BalancedSampling

import .Analysis: temporal_uncertainty


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

function quickplot(obs::Array, sim::Array, climate::Climate, label="", log=true; burn_in=1, limit=nothing, metric=Streamfall.mKGE)
    date = timesteps(climate)
    last_e = !isnothing(limit) ? limit : lastindex(obs)
    show_range = burn_in:last_e
    return quickplot(obs[show_range], sim[show_range], date[show_range], label, log; metric=metric)
end


function quickplot(obs::Array, sim::Array, xticklabels::Array, label="Modeled", log=true; metric=Streamfall.mKGE)
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

    combined = plot(fig, qqfig, size=(800, 400), left_margin=10mm, layout=(1,2))

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

https://discourse.julialang.org/t/symmetrical-log-plot/45709/3
"""
function symlog(y)
    y = replace(y, 0=>1e-06)  # small constant to avoid 0s
    # C = ceil(minimum(log10.(abs.(y))))
    return sign.(y) .* (log10.(abs.(y))) # / (10.0^C))
end


"""
    temporal_cross_section(dates, obs, sim; ylabel=nothing, func::Function=Streamfall.ME, period::Function=month)

Provides indication of predictive uncertainty across time, grouped by `period`.

Notes:
Assumes daily data.
Filters out leap days.

# Arguments
- dates : Date of each observation
- obs : observed data
- sim : modeled results
- ylabel : Optional replacement ylabel. Uses name of `func` if not provided.
- `func::Function` : Function to apply to each month-day grouping
- `period::Function` : Method from `Dates` package to group (defaults to `month`)
"""
function temporal_cross_section(dates, obs, sim; 
                                title="", ylabel=nothing, label=nothing, 
                                func::Function=Streamfall.ME, period::Function=monthday,
                                show_extremes::Bool=false, kwargs...)
    if isnothing(ylabel)
        ylabel = nameof(func)
    end

    if isnothing(label)
        label = ylabel
    end

    arg_keys = keys(kwargs)
    format_func = y -> y
    logscale = [:log, :log10]
    tmp = nothing
    if :yscale in arg_keys || :yaxis in arg_keys
        tmp = (:yscale in arg_keys) ? kwargs[:yscale] : kwargs[:yaxis]

        if tmp in logscale
            orig_obs = collect(obs)
            orig_sim = collect(sim)
            obs = symlog(obs)
            sim = symlog(sim)

            # Format function for y-axis tick labels (e.g., 10^x)
            format_func = y -> (y != 0) ? L"%$(round(sign(y) * 10))^{%$(round(abs(y), digits=1))}" : L"0"

            x_section, lower, upper, min_section, max_section, whisker_range, cv_r, std_error = temporal_uncertainty(dates, obs, sim, period, func)
            orig_x_section, _, _, _, _, _, _, orig_std_error = temporal_uncertainty(dates, orig_obs, orig_sim, period, func) 
        end
    else
        x_section, lower, upper, min_section, max_section, whisker_range, cv_r, std_error = temporal_uncertainty(dates, obs, sim, period, func)
    end

    sp = sort(unique(period.(dates)))
    deleteat!(sp, findall(x -> x == (2,29), sp))
    xlabels = join.(sp, "-")
    # whisker_range = upper .- lower

    if tmp in logscale
        # Remove keys 
        kwargs = Dict(kwargs)
        delete!(kwargs, :yscale)
        delete!(kwargs, :yaxis)

        # Display values using original data instead of log-transformed data
        m_ind = round(median(orig_x_section), digits=2)
        std_error = round(orig_std_error, digits=2)
    else
        m_ind = round(median(x_section), digits=2)
        std_error = round(std_error, digits=2)
    end 

    fig = plot(xlabels, x_section,
               label="$(label)\n[Median: $(m_ind) | SE: $(std_error)]",
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
               kwargs...)  # size=(1000, 350)

    plot!(fig, xlabels, lower, fillrange=upper, color="lightblue", alpha=0.5, label="", linealpha=0)

    if show_extremes
        scatter!(fig, xlabels, min_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
        scatter!(fig, xlabels, max_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0; kwargs...)
    end

    # plot!(fig, min_section, label="lowest", color="orange", alpha=0.5)
    # plot!(fig, max_section, label="highest", color="orange", alpha=0.5)

    # plot!(fig, lower, label="min q0.05", color="green", alpha=0.5)
    # plot!(fig, upper, label="max q0.95", color="green", alpha=0.5)

    return fig
end
