using Plots, StatsPlots
using Plots.Measures
using DataFrames, Dates, Statistics, Distributions
import StatsBase: ecdf
import Bootstrap: bootstrap, BalancedSampling


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
    date = timesteps(climate)
    last_e = !isnothing(limit) ? limit : lastindex(obs)
    show_range = burn_in:last_e
    return quickplot(obs[show_range], node.outflow[show_range], date[show_range], label, log; metric=metric)
end


function quickplot(obs::Array, sim::Array, xticklabels::Array, label="Modeled", log=true; metric=Streamfall.mKGE)
    @assert length(xticklabels) == length(obs) || "x-axis tick label length and observed lengths do not match!"
    @assert length(xticklabels) == length(sim) || "x-axis tick label length and simulated lengths do not match!"

    score = round(metric(obs, sim), digits=4)
    metric_name = String(Symbol(metric))

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
        # Add small constant in case of 0-flow
        obs = obs .+ 1e-2
        sim = sim .+ 1e-2
    end

    qqfig = qqplot(obs, sim, legend=false, markerstrokewidth=0, alpha=0.7, xlabel="Observed", ylabel="Modeled")

    # tick limits explicitly specified as workaround for known bug:
    # https://discourse.julialang.org/t/plotting-log-log-plot-how-do-i-resolve-a-no-strict-ticks-warning/67510
    if log
        xaxis!(qqfig, :log)
        # xlims!(minimum(sim)-10, maximum(sim)+10)
        yaxis!(qqfig, :log)
        # ylims!(minimum(sim)-10, maximum(sim)+10)
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
function temporal_cross_section(dates, obs, sim; title="", ylabel=nothing, func::Function=Streamfall.ME, period::Function=monthday)
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
        min_section[i] = Q1 - IQR_lim
        max_section[i] = Q3 + IQR_lim
        x_section[i] = Q2
    end

    if isnothing(ylabel)
        ylabel = nameof(func)
    end

    xlabels = join.(sp, "-")
    mean_ind = round(mean(x_section), digits=2)
    whisker_range = round(mean(max_section .- x_section), digits=2)
    fig = plot(xlabels, x_section,
               ribbon=(x_section .- min_section, max_section .- x_section),
               label="$(ylabel): [Mean: $(mean_ind) (± $(whisker_range))]",
               xlabel=nameof(period),
               ylabel=ylabel,
               legend=:bottomleft,
               legendfont=Plots.font(12),
               fg_legend=:transparent,
               bg_legend=:transparent,
               left_margin=5mm,
               bottom_margin=5mm,
               title=title,
               # yaxis=:log,
               size=(850, 400))

    return fig
end
