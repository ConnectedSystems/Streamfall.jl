using Statistics
using StatsBase


"""The Nash-Sutcliffe Efficiency score"""
NSE(obs, sim) = 1.0 - sum((obs .- sim).^2) / sum((obs .- mean(obs)).^2)


"""Normalized Nash-Sutcliffe Efficiency score (bounded between 0 and 1).

# References
1. Nossent, J., Bauwens, W., 2012.
    Application of a normalized Nash-Sutcliffe efficiency to improve the accuracy of the Sobol’ sensitivity analysis of a hydrological model.
    EGU General Assembly Conference Abstracts 237.
"""
NNSE(obs, sim) = 1.0 / (2.0 - NSE(obs, sim))


"""Root Mean Square Error"""
RMSE(obs, sim) = (sum((sim .- obs).^2)/length(sim))^0.5


"""Coefficient of determination (R^2)

Aliases `NSE()`
"""
function R2(obs, sim)::Float64
    return NSE(obs, sim)
end


"""Determine adjusted R^2

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
- `p::Int` : number of explanatory variables
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

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results

# References
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

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function BKGE(obs, sim)::Float64
    kge = KGE(obs, sim)
    return kge / (2 - kge)
end


"""Normalized KGE between 0 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function NKGE(obs, sim)::Float64
    return 1 / (2 - KGE(obs, sim))
end


"""Calculate the modified KGE metric (2012).

Also known as KGE prime (KGE').

# Arguments
- `obs::Vector`: observations
- `sim::Vector` : modeled results

# References
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


"""Bounded modified KGE between -1 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function BmKGE(obs, sim)::Float64
    mkge = mKGE(obs, sim)
    return mkge / (2 - mkge)
end


"""Normalized modified KGE between 0 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function NmKGE(obs, sim)::Float64
    return 1 / (2 - mKGE(obs, sim))
end


"""Calculate the non-parametric Kling-Gupta Efficiency (KGE) metric.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results

# References
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


"""Bounded non-parametric KGE between -1 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function BnpKGE(obs, sim)::Float64
    npkge = npKGE(obs, sim)
    return npkge / (2 - npkge)
end


"""Normalized non-parametric KGE between 0 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function NnpKGE(obs::Array, sim::Array)::Float64
    return 1 / (2 - npKGE(obs, sim))
end


"""Liu Mean Efficiency metric.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results

# References
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


function naive_split_metric(obs::Vector, sim::Vector, n_members::Int, metric::Function=NNSE)
    obs_chunks = collect(Iterators.partition(obs, n_members))
    sim_chunks = collect(Iterators.partition(sim, n_members))
    scores = Array{Float64, 1}(undef, length(obs_chunks))

    for (idx, h_chunk) in enumerate(obs_chunks)
        scores[idx] = metric(h_chunk, sim_chunks[idx])
    end

    return scores
end


"""Naive approach to split metrics.

Split metrics are a meta-objective optimization approach which segments data
into subperiods. The objective function is calculated for each subperiod and
then recombined. The approach addresses the lack of consideration of dry years 
with least-squares.

In Fowler et al., [1] the subperiod is one year. This method is "naive" in 
that the time series is partitioned into `N` chunks of `n_members`.
Therefore, leap years or partial years are not considered.

# Arguments
- `obs::Vector` : Historic observations to compare against
- `sim::Vector` : Modeled time series
- `n_members::Int` : number of members per chunk, defaults to 365
- `metric::Function` : Objective function to apply, defaults to NNSE
- `comb_method::Function` : Recombination method, defaults to `mean`

# References
1. Fowler, K., Peel, M., Western, A., Zhang, L., 2018.
    Improved Rainfall-Runoff Calibration for Drying Climate: Choice of Objective Function.
    Water Resources Research 54, 3392–3408.
    https://doi.org/10.1029/2017WR022466
"""
function naive_split_metric(obs, sim; n_members::Int=365, metric::Function=NNSE, comb_method::Function=mean)
    scores = naive_split_metric(obs, sim, n_members, metric)
    return comb_method(scores)
end