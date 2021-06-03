using Statistics
using StatsBase


NSE(obs, sim) = 1.0 - sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2)

NNSE(obs, sim) = 1.0 / (2.0 - NSE(obs, sim))

RMSE(obs, sim) = (sum((sim .- obs).^2)/length(sim))^0.5


"""Coefficient of determination (R^2)"""
function R2(obs, sim)::Float64
    return NSE(obs, sim)
end


"""Determine adjusted R^2

Parameters
----------
obs : observations
modeled : modeled results
p : number of explanatory variables
"""
function ADJ_R2(obs, sim, p::Int64)::Float64
    n = length(obs)
    adj_r2 = 1 - (1 - R2(obs, sim)) * ((n - 1) / (n - p - 1))

    return adj_r2
end


"""Calculate the 2009 Kling-Gupta Efficiency (KGE) metric.

A KGE score of 1 means perfect fit.
A score < -0.41 indicates that the mean of observations
provides better estimates (see Knoben et al., 2019).

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
function KGE(obs, sim)::Float64
    r = Statistics.cor(obs, sim)
    if isnan(r)
        r = 0.0
    end

    α = std(sim) / std(obs)
    # β = (mean(sim) - mean(obs)) / std(obs)
    β = mean(sim) / mean(obs)

    kge = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)

    return kge
end


"""Bounded KGE, bounded between -1 and 1.
"""
function BKGE(obs, sim)::Float64
    kge = KGE(obs, sim)
    return kge / (2 - kge)
end


"""Normalized KGE between 0 and 1"""
function NKGE(obs, sim)::Float64
    return 1 / (2 - KGE(obs, sim))
end


"""Calculate the modified KGE metric (2012).

Also known as KGE prime (KGE').

1. Kling, H., Fuchs, M., Paulin, M., 2012.
    Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios.
    Journal of Hydrology 424–425, 264–277.
    https://doi.org/10.1016/j.jhydrol.2012.01.011
"""
function mKGE(obs, sim)::Float64
    # Timing
    r = Statistics.cor(obs, sim)
    if isnan(r)
        r = 0.0
    end

    # Variability
    cv_s = StatsBase.variation(sim)
    if isnan(cv_s)
        cv_s = 0.0
    end
    γ = cv_s / StatsBase.variation(obs)

    # Magnitude
    β = mean(sim) / mean(obs)

    mod_kge = 1 - sqrt((r - 1)^2 + (β - 1)^2 + (γ - 1)^2)

    return mod_kge
end


"""Bounded modified KGE between -1 and 1."""
function BmKGE(obs, sim)::Float64
    mkge = mKGE(obs, sim)
    return mkge / (2 - mkge)
end


"""Normalized modified KGE between 0 and 1."""
function NmKGE(obs, sim)::Float64
    return 1 / (2 - mKGE(obs, sim))
end


"""Calculate the non-parametric Kling-Gupta Efficiency (KGE) metric.

References
----------
1. Pool, S., Vis, M., Seibert, J., 2018.
    Evaluating model performance: towards a non-parametric variant of the Kling-Gupta efficiency.
    Hydrological Sciences Journal 63, 1941–1953.
    https://doi.org/10.1080/02626667.2018.1552002

"""
function npKGE(obs, sim)::Float64
    r = StatsBase.corspearman(obs, sim)
    if isnan(r)
        r = 0.0
    end

    # flow duration curves
    μ_s = mean(sim)
    if μ_s == 0.0
        fdc_sim = repeat([0.0], length(sim))
    else
        x = length(sim) * μ_s
        fdc_sim = sort(sim / x)
    end

    μ_o = mean(obs)
    x = length(obs) * μ_o
    fdc_obs = sort(obs / x)

    α = 1 - 0.5 * sum(abs.(fdc_sim - fdc_obs))
    if μ_o == 0.0
        β = 0.0
    else
        β = μ_s / μ_o
    end

    kge = 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)

    return kge
end


"""Bounded non-parametric KGE between -1 and 1."""
function BnpKGE(obs, sim)::Float64
    npkge = npKGE(obs, sim)
    return npkge / (2 - npkge)
end


"""Normalized non-parametric KGE between 0 and 1."""
function NnpKGE(obs::Array, sim::Array)::Float64
    return 1 / (2 - npKGE(obs, sim))
end


"""Liu Mean Efficiency metric.

References
----------
1. Liu, D., 2020.
    A rational performance criterion for hydrological model.
    Journal of Hydrology 590, 125488.
    https://doi.org/10.1016/j.jhydrol.2020.125488

"""
function LME(obs, sim)::Float64
    μ_o = mean(obs)
    μ_s = mean(sim)
    β = (μ_s / μ_o)

    r = Statistics.cor(obs, sim)
    σ_s = std(sim)
    σ_o = std(obs)
    k_1 = r * (σ_s / σ_o)

    LME = 1 - sqrt((k_1 - 1)^2 + (β - 1)^2)

    return LME
end