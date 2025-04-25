module MakieExt

using Statistics
using DataFrames
using Dates
using LaTeXStrings

using Streamfall
import Streamfall: StreamfallNetwork
import Streamfall.Analysis: TemporalCrossSection

using Makie


function Streamfall.Viz.quickplot(node::NetworkNode)
    fig = lines(node.outflow)
    return fig
end

function Streamfall.Viz.quickplot(node::NetworkNode, climate::Climate)
    date = timesteps(climate)

    @assert length(date) == length(node.outflow) || "Date length and result lengths do not match!"

    f, ax, sp = lines(date, node.outflow)

    return f
end
function Streamfall.Viz.quickplot(node::DamNode, climate::Climate)
    date = timesteps(climate)

    @assert length(date) == length(node.level) || "Date length and result lengths do not match!"

    f, ax, sp = lines(date, node.level)

    return f
end

function Streamfall.Viz.quickplot(
    obs::Vector, node::NetworkNode, climate::Climate;
    label="Modeled", log=false, burn_in=1, limit=nothing, metric=Streamfall.mKGE
)
    return quickplot(obs, node.outflow, climate; label, log, burn_in=burn_in, limit=limit, metric=metric)
end
function Streamfall.Viz.quickplot(
    obs::Vector, node::DamNode, climate::Climate;
    label="Modeled", log=false, burn_in=1, limit=nothing, metric=Streamfall.mKGE
)
    return quickplot(obs, node.level, climate; label, log, burn_in=burn_in, limit=limit, metric=metric)
end

function Streamfall.Viz.quickplot(
    obs::DataFrame, node::NetworkNode, climate::Climate;
    label="", log=false, burn_in=1, limit=nothing, metric=Streamfall.mKGE
)
    return Streamfall.Viz.quickplot(obs[:, node.name], node.outflow, climate; label, log, burn_in, limit, metric)
end
function Streamfall.Viz.quickplot(
    obs::DataFrame, node::DamNode, climate::Climate;
    label="", log=false, burn_in=1, limit=nothing, metric=Streamfall.mKGE
)
    return Streamfall.Viz.quickplot(obs[:, node.name], node.level, climate; label, log, burn_in, limit, metric)
end

function Streamfall.Viz.quickplot(
    obs::Vector, sim::Vector, climate::Climate;
    label="Modeled", log=false, burn_in=1, limit=nothing, metric=Streamfall.mKGE
)
    date = timesteps(climate)
    last_e = !isnothing(limit) ? limit : lastindex(obs)
    show_range = burn_in:last_e
    return Streamfall.Viz.quickplot(obs[show_range], sim[show_range], date[show_range], label, log; metric=metric)
end
function Streamfall.Viz.quickplot(obs::Vector, sim::Vector, xticklabels::Vector, label="Modeled", log=false; metric=Streamfall.mKGE)
    @assert length(xticklabels) == length(obs) || "x-axis tick label length and observed lengths do not match!"
    @assert length(xticklabels) == length(sim) || "x-axis tick label length and simulated lengths do not match!"

    score = round(metric(obs, sim), digits=4)
    metric_name = String(Symbol(metric))

    if log
        # Add small constant in case of 0-flow
        obs = copy(obs)
        sim = copy(sim)
        obs[obs.==0.0] .+= 1e-4
        sim[sim.==0.0] .+= 1e-4
    end

    scale = log == false ? identity : log10

    f = Figure(size=(850, 400))
    flow_ax = Axis(f[1, 1]; xlabel="Date", ylabel="Streamflow", yscale=scale)
    qq_ax = Axis(f[1, 2]; xlabel="Observed", ylabel="Modeled", xscale=scale, yscale=scale)

    label = "$(label) ($(metric_name): $(score))"
    lines!(flow_ax, xticklabels, obs, label="Observed")
    lines!(flow_ax, xticklabels, sim; label=label, alpha=0.5)

    qqplot!(
        qq_ax, obs[burn_in:end], sim[burn_in:end];
        qqline=:identity, strokewidth=0.01, strokecolor=(:blue, 0.01), # legend=false, strokewidth=0.03, strokealpha=0.1,
        markercolor=(:blue, 0.02)
    )
    leg = Legend(f[2, 1:2], flow_ax; orientation=:horizontal)

    return f
end

"""
    plot_residuals(obs::AbstractVector, sim::AbstractVector; xlabel="", ylabel="", title="")

Plot residual between two sequences.

# Arguments
- x : x-axis data
- y : y-axis data
- xlabel : x-axis label
- ylabel : y-axis label
- title : title text
"""
function Streamfall.Viz.plot_residuals(x::AbstractVector, y::AbstractVector; xlabel="", ylabel="", title="")
    # 1:1 Plot
    f, ax, sp = scatter(x, y, strokewidth=0, alpha=0.2)

    # Calculate means
    x_mean = mean(x)
    y_mean = mean(y)

    # Calculate slope (β)
    β = sum((x .- x_mean) .* (y .- y_mean)) / sum((x .- x_mean) .^ 2)

    # Calculate intercept (α)
    α = y_mean - β * x_mean

    # Approximate x values for the line
    x_line = range(minimum(x), maximum(x), length=length(x))

    # Calculate corresponding y values
    y_line = α .+ β .* x_line

    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.title = title

    lines!(ax, y_line; color=(:black, 0.3))

    return f
end

"""
    temporal_cross_section(
        dates, obs;
        title="",
        ylabel="Mean",
        label=nothing,
        period=monthday,
        kwargs...
    )

Provides indication of temporal variation and uncertainty across time, grouped by `period`.

Notes:
Assumes daily data.
Filters out leap days.

# Arguments
- `dates` : Date of each observation
- `obs` : observed data
- `title` : Plot title
- `ylabel` : Optional replacement ylabel. Uses name of `func` if not provided.
- `period` : Method from `Dates` package to group (defaults to `monthday`)
- `kwargs` : Additional plotting keyword arguments
"""
function Streamfall.Viz.temporal_cross_section(
    dates, obs;
    title="", ylabel="Mean", label=nothing,
    period::Function=monthday,
    kwargs...
)
    if isnothing(label)
        label = ylabel
    end

    arg_keys = keys(kwargs)
    format_func = y -> y
    logscale = [:log, :log10]
    tmp = nothing

    xsect_res = TemporalCrossSection(dates, obs, period)

    # TODO: This is currently broken
    if :yscale in arg_keys || :yaxis in arg_keys
        tmp = (:yscale in arg_keys) ? kwargs[:yscale] : kwargs[:yaxis]

        if tmp in logscale
            log_obs = symlog(copy(obs))
            format_func = y -> (y != 0) ? L"%$(Int(round(sign(y)) * 10))^{%$(round(abs(y), digits=1))}" : L"0"
            log_xsect_res = TemporalCrossSection(dates, log_obs, period)
            target = log_xsect_res.cross_section
        else
            target = xsect_res.cross_section
        end
    else
        target = xsect_res.cross_section
    end

    x_section = target.mean
    lower_75 = target.lower_75
    upper_75 = target.upper_75
    lower_95 = target.lower_95
    upper_95 = target.upper_95

    sp = target.subperiod
    xlabels = join.(sp, "-")

    # Calculate statistics for legend
    xs_mean = xsect_res.cross_section.mean
    m_ind = round(mean(xs_mean), digits=2)
    sd_ind = round(std(xs_mean), digits=2)
    wr75_m_ind = round(xsect_res.mean_75, digits=2)
    wr75_sd_ind = round(xsect_res.std_75, digits=2)
    wr95_m_ind = round(xsect_res.mean_95, digits=2)
    wr95_sd_ind = round(xsect_res.std_95, digits=2)

    # Create figure
    fig = Figure(size=(800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel=string(nameof(period)),
        ylabel=ylabel,
        title=title
    )

    # Plot confidence intervals and median line
    x = 1:length(xlabels)

    # 95% CI
    b95 = band!(
        ax, x, lower_95, upper_95,
        color=(:lightblue, 0.3)
    )

    # 75% CI
    b75 = band!(
        ax, x, lower_75, upper_75,
        color=(:lightblue, 0.5)
    )

    # Median line
    med_line = lines!(
        ax, x, x_section,
        color=:black
    )

    # Set x-axis labels
    if period == monthday
        ax.xticks = (1:14:length(xlabels), xlabels[1:14:end])
    elseif period == yearmonth
        ax.xticks = (1:12:length(xlabels), xlabels[1:12:end])
    end

    ax.xticklabelrotation = π / 4

    # Apply log scale if specified
    if !isnothing(tmp) && (tmp in logscale)
        ax.yscale = log10
        ax.formatter = format_func
    end

    # Add legend at bottom
    Legend(
        fig[2, 1],
        [med_line, b75, b95],
        [
            "Mean of $(label) μ: $(m_ind), σ: $(sd_ind)",
            "CI₇₅ μ: $(wr75_m_ind), σ: $(wr75_sd_ind)",
            "CI₉₅ μ: $(wr95_m_ind), σ: $(wr95_sd_ind)"
        ],
        orientation=:horizontal,
        tellwidth=false,
        tellheight=true
    )

    # Adjust layout
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))

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
function Streamfall.Viz.temporal_cross_section(
    dates, obs, sim;
    title="", ylabel="Median Error", label=nothing, period::Function=monthday, kwargs...
)
    if isnothing(label)
        label = ylabel
    end

    arg_keys = keys(kwargs)
    format_func = y -> y
    logscale = [:log, :log10]
    tmp = nothing

    xsect_res = TemporalCrossSection(dates, obs, sim, period)
    target = xsect_res.cross_section

    # TODO: This is currently broken
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

    xsect_df = xsect_res.cross_section

    # Calculate statistics for legend
    med = xsect_df.median
    m_ind = round(mean(med), digits=2)
    sd_ind = round(std(med), digits=2)
    wr75_m_ind = round(xsect_res.mean_75, digits=2)
    wr75_sd_ind = round(xsect_res.std_75, digits=2)
    wr95_m_ind = round(xsect_res.mean_95, digits=2)
    wr95_sd_ind = round(xsect_res.std_95, digits=2)

    x = 1:length(xlabels)

    # Create figure
    fig = Figure(size=(800, 600))
    ax = Axis(
        fig[1, 1],
        xlabel=string(nameof(period)),
        ylabel=ylabel,
        title=title
    )

    # Plot confidence intervals and median line
    # 95% CI
    b95 = band!(
        ax, x, lower_95, upper_95,
        color=(:lightblue, 0.3),
        label="CI₉₅ μ: $(wr95_m_ind), σ: $(wr95_sd_ind)"
    )

    # 75% CI
    b75 = band!(
        ax, x, lower_75, upper_75,
        color=(:lightblue, 0.5),
        label="CI₇₅ μ: $(wr75_m_ind), σ: $(wr75_sd_ind)"
    )

    # Median line
    med_line = lines!(
        ax, x, x_section,
        color=:black,
        label="$(label) μ: $(m_ind), σ: $(sd_ind)"
    )

    # Set x-axis labels
    if period == monthday
        ax.xticks = (1:14:length(xlabels), xlabels[1:14:end])
    elseif period == yearmonth
        ax.xticks = (1:12:length(xlabels), xlabels[1:12:end])
    end

    ax.xticklabelrotation = π / 4

    # Apply log scale if specified
    if !isnothing(tmp) && (tmp in logscale)
        ax.yscale = log10
        ax.formatter = format_func
    end

    # Add legend at bottom
    Legend(
        fig[2, 1],
        [med_line, b75, b95],
        [
            "$(label) μ: $(m_ind), σ: $(sd_ind)",
            "CI₇₅ μ: $(wr75_m_ind), σ: $(wr75_sd_ind)",
            "CI₉₅ μ: $(wr95_m_ind), σ: $(wr95_sd_ind)"
        ],
        orientation=:horizontal,
        tellwidth=false,
        tellheight=true,
        framevisible=false
    )  # Equivalent to transparent background

    # Adjust layout
    rowsize!(fig.layout, 1, Relative(0.9))
    rowsize!(fig.layout, 2, Relative(0.1))

    # Process any additional keyword arguments
    for (key, value) in kwargs
        if key ∉ [:yscale, :yaxis]  # Skip already processed arguments
            setproperty!(ax, key, value)
        end
    end

    return fig
end

"""
    save_figure(f::Figure, fn::String)

Save a figure.
"""
function Streamfall.Viz.save_figure(f::Figure, fn::String)
    save(fn, f, update=false)
end
function Streamfall.Viz.save_figure!(fn::String)
    f = current_figure()
    Streamfall.Viz.save_figure(f, fn)
end

end