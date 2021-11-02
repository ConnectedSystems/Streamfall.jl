using Plots, StatsPlots
using Plots.Measures
using DataFrames, Dates, Statistics, Distributions
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
                                show_extremes::Bool=false)
    x_section, lower, upper, min_section, max_section, whisker_range, cv_r = temporal_uncertainty(dates, obs, sim, period, func)

    if isnothing(ylabel)
        ylabel = nameof(func)
    end

    if isnothing(label)
        label = ylabel
    end

    sp = sort(unique(period.(dates)))
    deleteat!(sp, findall(x -> x == (2,29), sp))
    xlabels = join.(sp, "-")
    m_ind = round(median(x_section), digits=2)
    whisker_range = upper .- lower
    r_ind = round(mean(whisker_range), digits=2)

    cv_r = round(cv_r, digits=2)

    fig = plot(xlabels, x_section,
               ribbon=(x_section .+ abs.(lower), upper .- x_section),
               label="$(label)\n[Median: $(m_ind) | Mean CIᵣ: $(r_ind) | CVᵣ: $(cv_r)]",
               xlabel=nameof(period),
               ylabel=ylabel,
               legend=:bottomleft,
               legendfont=Plots.font(10),
               fg_legend=:transparent,
               bg_legend=:transparent,
               left_margin=5mm,
               bottom_margin=5mm,
               title=title,
               size=(1000, 350))

    if show_extremes
        scatter!(xlabels, max_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0)
        scatter!(xlabels, min_section, label="", alpha=0.5, color="lightblue", markerstrokewidth=0)
    end

    return fig
end
