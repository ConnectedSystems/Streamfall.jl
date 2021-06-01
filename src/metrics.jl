using Statistics
using StatsBase


NSE(obs, modeled) = 1.0 - sum((obs .- modeled).^2) / sum((obs .- mean(obs)).^2)

NNSE(obs, modeled) = 1.0 / (2.0 - NSE(obs, modeled))

RMSE(obs, modeled) = (sum((modeled .- obs).^2)/length(modeled))^0.5


"""Coefficient of determination (R^2)"""
function R2(obs::Array, modeled::Array)::Float64
    ss_tot = sum((obs .- mean(obs)).^2)
    ss_res = sum((obs .- modeled).^2)

    local r2 = 1 - (ss_res / ss_tot)
    return r2
end


"""Determine adjusted R^2

Parameters
----------
obs : observations
modeled : modeled results
p : number of explanatory variables
"""
function ADJ_R2(obs::Array, modeled::Array, p::Int64)::Float64
    n = length(obs)
    adj_r2 = 1 - (1 - R2(obs, modeled)) * ((n - 1) / (n - p - 1))

    return adj_r2
end


"""Calculate the 2009 Kling-Gupta Efficiency (KGE) metric.

A KGE score of 1 means perfect fit.
Similar to NSE, a score < 0 indicates that a mean of observations 
provides better estimates (although this is debated, see [2]).

Note: Although similar, NSE and KGE cannot be directly compared.

References
----------
1. Gupta, H.V., Kling, H., Yilmaz, K.K., Martinez, G.F., 2009. 
    Decomposition of the mean squared error and NSE performance criteria: 
    Implications for improving hydrological modelling. 
    Journal of Hydrology 377, 80–91. 
    https://doi.org/10.1016/j.jhydrol.2009.08.003

2. Knoben, W.J.M., Freer, J.E., Woods, R.A., 2019. 
    Technical note: Inherent benchmark or not? Comparing Nash-Sutcliffe and Kling-Gupta efficiency scores (preprint). 
    Catchment hydrology/Modelling approaches. 
    https://doi.org/10.5194/hess-2019-327

"""
function KGE(obs::Array{Float64}, modeled::Array{Float64})::Float64
    r = Statistics.cor(obs, modeled)
    α = std(modeled) / std(obs)
    β = (mean(modeled) - mean(obs)) / std(obs)    

    kge = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)

    return kge
end


"""Normalized KGE, bounded between 0 and 1.
"""
function NKGE(obs::Array{Float64}, sim::Array{Float64}) 
    kge = KGE(obs, sim)
    return kge / (2 - kge)
end


"""Calculate the modified KGE metric (2012).

1. Kling, H., Fuchs, M., Paulin, M., 2012. 
    Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios. 
    Journal of Hydrology 424–425, 264–277. 
    https://doi.org/10.1016/j.jhydrol.2012.01.011
"""
function mKGE(obs::Array{Float64}, sim::Array{Float64})::Float64
    r = Statistics.cor(obs, sim)
    γ = StatsBase.variation(sim) / StatsBase.variation(obs)
    β = (mean(sim) - mean(obs)) / std(obs)    

    mod_kge = 1 - sqrt((r - 1)^2 + (β - 1)^2 + (γ - 1)^2)

    return mod_kge
end


"""Normalized modified KGE, bounded between 0 and 1."""
function NmKGE(obs::Array{Float64}, sim::Array{Float64}) 
    mkge = mKGE(obs, sim)
    return mkge / (2 - mkge)
end


"""Calculate the non-parametric Kling-Gupta Efficiency (KGE) metric.

References
----------
1. Pool, S., Vis, M., Seibert, J., 2018.
    Evaluating model performance: towards a non-parametric variant of the Kling-Gupta efficiency. 
    Hydrological Sciences Journal 63, 1941–1953. 
    https://doi.org/10.1080/02626667.2018.1552002

"""
function npKGE(obs::Array{Float64}, sim::Array{Float64})::Float64
    r = StatsBase.corspearman(obs, sim)

    x = length(sim) * mean(sim)
    fdc_sim = sort(sim / x)
    
    x = length(obs) * mean(obs)
    fdc_obs = sort(obs / x)

    α = 1 - 0.5 * sum(abs(fdc_sim - fdc_obs))
    β = (mean(sim) - mean(obs)) / std(obs)    

    kge = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)

    return kge
end


"""Normalized non-parametric KGE, bounded between 0 and 1."""
function NnpKGE(obs::Array{Float64}, sim::Array{Float64})::Float64
    npkge = npKGE(obs, sim)
    return npkge / (2 - npkge)
end