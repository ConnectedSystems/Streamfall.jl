using Plots, StatsPlots


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


function quickplot(obs, node::NetworkNode, climate::Climate, label=""; burn_in=1, limit=nothing)
    date = timesteps(climate)
    last_e = !isnothing(limit) ? limit : lastindex(obs)
    show_range = burn_in:last_e
    return quickplot(obs[show_range], node.outflow[show_range], date[show_range], label)
end


function quickplot(obs::Array, sim::Array, xticklabels::Array, label="")
    @assert length(xticklabels) == length(obs) || "x-axis tick label length and observed lengths do not match!"
    @assert length(xticklabels) == length(sim) || "x-axis tick label length and simulated lengths do not match!"
    
    fig = plot(xticklabels, obs, label="Historic")
    plot!(xticklabels, sim, label=label)
    qqfig = qqplot(obs, sim, legend=false, markerstrokewidth=0)

    combined = plot(fig, qqfig, size=(800, 400))

    return combined
end
