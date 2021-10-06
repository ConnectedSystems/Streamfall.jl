using Statistics
using StatsBase


"""
Bounds given metric between -1.0 and 1.0, where 1.0 is perfect fit.

Suitable for use with any metric that ranges from 1 to -∞.

# References
1. Mathevet, T., Michel, C., Andréassian, V., Perrin, C., 2006.
    A bounded version of the Nash-Sutcliffe criterion for better model
    assessment on large sets of basins.
    IAHS-AISH Publication 307, 211–219.
    https://iahs.info/uploads/dms/13614.21--211-219-41-MATHEVET.pdf

# Example
```julia
julia> import Streamfall: @bound, KGE
julia> @bound KGE([1,2], [3,2])
-0.35653767993482094
```
"""
macro bound(metric)
    tmp = :($metric)
    return :($tmp / (2.0 - $tmp))
end


"""
Normalizes given metric between 0.0 and 1.0, where 1.0 is perfect fit.

Suitable for use with any metric that ranges from 1 to -∞.

# References
1. Nossent, J., Bauwens, W., 2012.
    Application of a normalized Nash-Sutcliffe efficiency to improve the
    accuracy of the Sobol’ sensitivity analysis of a hydrological model.
    EGU General Assembly Conference Abstracts 237.

# Example
```julia
julia> import Streamfall: @normalize, KGE
julia> @normalize KGE([1,2], [3,2])
0.1111111111111111
```
"""
macro normalize(metric)
    return :(1.0 / (2.0 - $metric))
end


"""
Applies mean inverse approach to a metric.

Suitable for use with any metric that ranges from 1 to -∞.

If using with other macros such as `@normalize` or `@bound`,
these must come first.

# References
1. Garcia, F., Folton, N., Oudin, L., 2017.
    Which objective function to calibrate rainfall–runoff
        models for low-flow index simulations?
    Hydrological Sciences Journal 62, 1149–1166.
    https://doi.org/10.1080/02626667.2017.1308511

# Example
```julia
julia> import Streamfall: @normalize, @mean_inverse, KGE
julia> @normalize @mean_inverse KGE [1,2] [3,2] 1e-6
0.3193506006429825
```
"""
macro mean_inverse(metric, obs, sim, ϵ=1e-2)
    obj, o, s, ϵ = eval(metric), eval(obs), eval(sim), eval(ϵ)
    q = obj(o, s)
    q2 = obj(1.0 ./ (o .+ ϵ), 1.0 ./ (s .+ ϵ))
    return mean([q, q2])
end


"""
Applies split meta metric approach

If using with other macros such as `@normalize` or `@bound`,
these must come first.

# References
1. Fowler, K., Peel, M., Western, A., Zhang, L., 2018.
    Improved Rainfall-Runoff Calibration for Drying Climate:
    Choice of Objective Function.
    Water Resources Research 54, 3392–3408.
    https://doi.org/10.1029/2017WR022466

# Example
```julia
julia> using Statistics
julia> import Streamfall: @normalize, @split, KGE
julia> @normalize @split KGE repeat([1,2], 365) repeat([3,2], 365) 365 mean
0.3217309561946589
```
"""
macro split(metric, obs, sim, n, agg_func=mean)
    obj, o, s, t, func = eval(metric), eval(obs), eval(sim), eval(n), eval(agg_func)
    return naive_split_metric(o, s; n_members=t, metric=obj, comb_method=func)
end


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


"""Adjusted R^2

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


"""
Mean Absolute Error
"""
MAE(obs, sim) = mean(abs.(sim .- obs))


"""
Mean Error
"""
ME(obs, sim) = mean(sim .- obs)


"""
    PBIAS(obs::Vector, sim::Vector)::Float64

Percent bias between `sim` and `obs`

Model performance for streamflow can be determined to be
satisfactory if the Nash-Sutcliffe Efficiency (NSE)
score > 0.5, the RMSE standard deviation ratio (RSR) < 0.7
and percent bias (PBIAS) is +/- 25% (see [1]).

# References
1. Moriasi, D.N., Arnold, J.G., Liew, M.W.V., Bingner, R.L.,
    Harmel, R.D., Veith, T.L., 2007.
    Model Evaluation Guidelines for Systematic Quantification
    of Accuracy in Watershed Simulations.
    Transactions of the ASABE 50, 885–900.
    https://doi.org/10.13031/2013.23153
"""
PBIAS(obs, sim) = (sum(obs .- sim) * 100) / sum(obs)


"""
    RSR(obs::Vector, sim::Vector)::Float64

The RMSE-observations standard deviation ratio (RSR).

Varies between 0 and a large positive value, where 0
indicates an RMSE value of 0.

# References
1. Moriasi, D.N., Arnold, J.G., Liew, M.W.V., Bingner, R.L.,
    Harmel, R.D., Veith, T.L., 2007.
    Model Evaluation Guidelines for Systematic Quantification
    of Accuracy in Watershed Simulations.
    Transactions of the ASABE 50, 885–900.
    https://doi.org/10.13031/2013.23153
"""
function RSR(obs, sim)::Float64
    rmse = RMSE(obs, sim)
    σ_obs = std(obs)
    rsr = rmse / σ_obs
    return rsr
end


"""
    EV(obs, sim)

Explained Variance.

Represents the amount of variation in the observations which the predictions
are able to explain.
"""
function EV(obs, sim)
    return 1.0 - (var(obs - sim) / obs)
end


"""
Relative Skill Score.

Provides an indication of model performance relative to a known benchmark score.

Suitable for use with least-squares approaches that provide skill scores ranging from 1 to -∞.

# Arguments
- `Sb` : Benchmark score
- `Sm` : Model score

# References
1. Knoben, W.J.M., Freer, J.E., Woods, R.A., 2019.
    Technical note: Inherent benchmark or not?
        Comparing Nash-Sutcliffe and Kling-Gupta efficiency scores (preprint).
    Catchment hydrology/Modelling approaches.
    https://doi.org/10.5194/hess-2019-327
"""
function relative_skill_score(Sb::Float64, Sm::Float64)
    return (Sm - Sb) / (1.0 - Sb)
end


"""
    NSE_logbias(obs, sim; metric::Function=NSE, bias_threshold::Float64=5.0, shape::Float64=2.5)

The NSE_logbias meta-metric provides a weighted combination of a least-squares approach and a
logarithmic function of bias. The metric penalizes predictions with an overall bias above a
threshold (defined as 5% in [1]).

It is also referred to as the Viney F score.

# Extended help
The penalty applied is non-symmetrical (or multiplicatively symmetrical) in that
predictions that are _double_ the observed are penalized identically to predictions that are
_half_ the observed volume.

# Arguments
- `obs::Vector` : Historic observations to compare against
- `sim::Vector` : Modeled time series
- `metric::Function` : least-squares method to use, defaults to NSE
- `bias_threshold::Float64` : Bias threshold after which the score given by `metric` is penalized, defaults to 5 (%)
- `shape::Float64` : Exponent value controlling shape of penalization (see Figure 2 in [1]).

# References
1. Viney, N. R., Perraud, J., Vaze, J., Chiew, F.H.S., Post, D.A., Yang, A. 2009
    The usefulness of bias constraints in model calibration for regionalisation
        to ungauged catchments
    18th World IMACS / MODSIM Congress, Cairns, Australia, 13 - 17 July 2009
    Available at:
    https://www.researchgate.net/publication/294697092_The_usefulness_of_bias_constraints_in_model_calibration_for_regionalisation_to_ungauged_catchments

2. Teng, J., Potter, N.J., Chiew, F.H.S., Zhang, L., Wang, B., Vaze, J., Evans, J.P., 2015.
    How does bias correction of regional climate model precipitation affect modelled runoff?
    Hydrology and Earth System Sciences 19, 711–728.
    https://doi.org/10.5194/hess-19-711-2015
"""
function NSE_logbias(obs, sim; metric=NSE, bias_threshold=5.0, shape=2.5)
    E_ns = metric(obs, sim)
    B = PBIAS(obs, sim) / 100.0
    return E_ns - bias_threshold * abs(log(1 + B))^shape
end


"""
    KGE(obs::Vector, sim::Vector; scaling::Tuple=nothing)::Float64

Calculate the 2009 Kling-Gupta Efficiency (KGE) metric.

Decomposes NSE into correlation (`r`), relative variability (`α`), and bias (`β`) terms.

A KGE score of 1 means perfect fit.
A score < -0.41 indicates that the mean of observations
provides better estimates (see Knoben et al., [2]).

The `scaling` argument expects a three-valued tuple
which scales `r`, `α` and `β` factors respectively.
If not specified, defaults to `1`.

Note: Although similar, NSE and KGE cannot be directly compared.

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

3. Mizukami, N., Rakovec, O., Newman, A.J., Clark, M.P., Wood, A.W.,
    Gupta, H.V., Kumar, R., 2019.
    On the choice of calibration metrics for “high-flow”
        estimation using hydrologic models.
    Hydrology and Earth System Sciences 23, 2601–2614.
    https://doi.org/10.5194/hess-23-2601-2019
"""
function KGE(obs, sim; scaling=nothing)::Float64
    if isnothing(scaling)
        scaling = (1, 1, 1)
    end

    # Correlation
    r = Statistics.cor(obs, sim)
    if isnan(r)
        r = 0.0
    end

    # relative variance
    α = std(sim) / std(obs)

    # bias
    β = mean(sim) / mean(obs)

    rs = scaling[1]
    as = scaling[2]
    bs = scaling[3]

    kge = 1 - sqrt(rs*(r - 1)^2 + as*(α - 1)^2 + bs*(β - 1)^2)

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
function NKGE(obs, sim; scaling=nothing)::Float64
    return 1 / (2 - KGE(obs, sim; scaling=scaling))
end


"""Calculate the modified KGE metric (2012).

Also known as KGE prime (KGE').

# Extended help

It is not recommended to apply KGE' with log-transformed flows (see [2]).
Numerical instabilities arise as flow approaches values close to zero.
This is possible under extreme dry conditions or by chance when sub-sampling.

In cases where observations are constant or otherwise displays zero variance or zero
mean flow, this implementation applies a simple logistic function (ℯ⁻ˣ) to gain an 
indication of simulated data's distance to zero.

This is to:
- avoid NaNs influencing subsequent calculations 
- allow use with split methods which may partition streamflows into periods of 0 flows.

# Arguments
- `obs::Vector`: observations
- `sim::Vector` : modeled results
- `scaling::Tuple` : scaling factors in order of timing (r), magnitude (β), variability (γ).
                     Defaults to (1,1,1).

# References
1. Kling, H., Fuchs, M., Paulin, M., 2012.
    Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios.
    Journal of Hydrology 424–425, 264–277.
    https://doi.org/10.1016/j.jhydrol.2012.01.011

2. Santos, L., Thirel, G., Perrin, C., 2018.
    Technical note: Pitfalls in using log-transformed flows within the KGE criterion.
    Hydrology and Earth System Sciences 22, 4583–4591.
    https://doi.org/10.5194/hess-22-4583-2018
"""
function mKGE(obs, sim; scaling=nothing)::Float64
    if isnothing(scaling)
        scaling = (1,1,1)
    end

    # Timing (Pearson's correlation)
    r = Statistics.cor(obs, sim)
    if isnan(r)
        # can happen if obs or sim is a constant (std of 0)
        r = all(obs .== sim) ? 1.0 : 0.0
    end

    # Variability
    cv_s = StatsBase.variation(sim)
    if isnan(cv_s)
        cv_s = 0.0
    end

    cv_o = StatsBase.variation(obs)
    if isnan(cv_o)
        cv_o = 0.0
    end

    if cv_o == 0.0 && cv_s == 0.0
        γ = 1.0
    elseif cv_o == 0.0
        # Use logistic function to indicate distance from 0
        γ = 1.0 / exp(cv_s)
    else
        γ = cv_s / cv_o
    end

    # Magnitude
    μ_o = mean(obs)
    μ_s = mean(sim)
    if μ_o == 0.0
        # use logistic function to indicate distance from 0 (μ_o)
        β = 1.0 / exp(μ_s)
    else
        β = μ_s / μ_o
    end

    rs = scaling[1]
    βs = scaling[2]
    γs = scaling[3]

    mod_kge = 1.0 - sqrt(rs*(r - 1)^2 + βs*(β - 1)^2 + γs*(γ - 1)^2)

    return mod_kge
end


"""Bounded modified KGE between -1 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function BmKGE(obs, sim; scaling=nothing)::Float64
    mkge = mKGE(obs, sim; scaling=scaling)
    return mkge / (2 - mkge)
end


"""Normalized modified KGE between 0 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function NmKGE(obs, sim; scaling=nothing)::Float64
    return 1 / (2 - mKGE(obs, sim; scaling=scaling))
end


"""
Mean Inverse NmKGE

Said to produce better fits for low-flow indices
compared to mKGE (see [1]).

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
- `scaling::Tuple` : scaling factors for r, α, and β (defaults to 1.0)
- `ϵ::Float64` : small constant to use with inverse flow to allow consideration of periods with no flow. Defaults to 1e-2.

# References
1. Garcia, F., Folton, N., Oudin, L., 2017.
    Which objective function to calibrate rainfall–runoff
        models for low-flow index simulations?
    Hydrological Sciences Journal 62, 1149–1166.
    https://doi.org/10.1080/02626667.2017.1308511
"""
mean_NmKGE(obs, sim; scaling=nothing, ϵ=1e-2) = mean([Streamfall.NmKGE(obs, sim; scaling=scaling), Streamfall.NmKGE(1.0 ./ (obs .+ ϵ), 1.0 ./ (sim .+ ϵ); scaling=scaling)])


"""Calculate the non-parametric Kling-Gupta Efficiency (KGE) metric.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled
- `scaling::Tuple` : scaling factors for timing (s), variability (α), magnitude (β)

# References
1. Pool, S., Vis, M., Seibert, J., 2018.
    Evaluating model performance: towards a non-parametric variant of the Kling-Gupta efficiency.
    Hydrological Sciences Journal 63, 1941–1953.
    https://doi.org/10.1080/02626667.2018.1552002

"""
function npKGE(obs, sim; scaling=nothing)::Float64
    if isnothing(scaling)
        scaling = (1,1,1)
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
    if μ_o == 0.0
        fdc_obs = repeat([0.0], length(obs))
    else
        x = length(obs) * μ_o
        fdc_obs = sort(obs / x)
    end

    α = 1 - 0.5 * sum(abs.(fdc_sim - fdc_obs))

    # Magnitude
    if μ_o == 0.0
        β = μ_s == 0.0 ? 1.0 : μ_s
    else
        β = μ_s / μ_o
    end

    # Timing and flow dynamics
    r = StatsBase.corspearman(fdc_obs, fdc_sim)
    if isnan(r)
        r = 1.0  # can occur if identical sequences are used (e.g., 0 flows)
    end

    rs = scaling[1]
    αs = scaling[2]
    βs = scaling[3]

    kge = 1 - sqrt(rs*(r - 1)^2 + αs*(α - 1)^2 + βs*(β - 1)^2)

    return kge
end


"""Bounded non-parametric KGE between -1 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function BnpKGE(obs, sim; scaling=nothing)::Float64
    npkge = npKGE(obs, sim; scaling=scaling)
    return npkge / (2 - npkge)
end


"""Normalized non-parametric KGE between 0 and 1.

# Arguments
- `obs::Vector` : observations
- `sim::Vector` : modeled results
"""
function NnpKGE(obs, sim; scaling=nothing)::Float64
    return 1 / (2 - npKGE(obs, sim; scaling=scaling))
end


"""Liu Mean Efficiency metric (LME).

Reformulation of the KGE metric said to be advantageous for capturing extreme
flow events.

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
    α = (σ_s / σ_o)
    k_1 = r * α

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

Split metrics are a meta-objective optimization approach which "splits" data
into subperiods. The objective function is calculated for each subperiod and
then recombined. The approach addresses the lack of consideration of dry years
with least-squares.

In Fowler et al., [1] the subperiod is one year. The implementation offered here
is "naive" in that the data is partitioned into `N` chunks of `n_members` and
does not consider date/time.

# Arguments
- `obs::Vector` : Historic observations to compare against
- `sim::Vector` : Modeled time series
- `n_members::Int` : number of members per chunk (i.e., sub-samples), defaults to 365
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


"""
    inverse_metric(obs, sim; metric, comb_method::Function=mean)

A meta-objective function which combines the performance of the
given metric as applied to the discharge and the inverse of the
discharge.

By default, the combination method is to take the mean.

# Arguments
- obs : observed
- sim : modeled results
- `metric::Function` : objective function
- `comb_method::Function` : mean
- ϵ : offset value to use (enables use with zero-flow time steps), defaults to 1e-2


# References
1. Garcia, F., Folton, N., Oudin, L., 2017.
    Which objective function to calibrate rainfall–runoff models
        for low-flow index simulations?
    Hydrological Sciences Journal 62, 1149–1166.
    https://doi.org/10.1080/02626667.2017.1308511
"""
function inverse_metric(obs, sim, metric::Function; comb_method::Function=mean, ϵ=1e-2)
    return comb_method([metric(obs, sim), metric(1.0 ./ (obs + ϵ), 1.0 ./ (sim + ϵ))])
end


"""
Allows comparison of any model compared against a pre-defined benchmark, assuming both scores were obtained with
the same objective function.

Positive values indicate a model is better than the benchmark, and negative values indicate a model performs worse.

# Extended help
It is noted in Knoben et al., [1] that the skill score should always be contextualized with the original benchmark
value. Interpreting skill scores by themselves may become difficult if the benchmark score is already quite high.
A small improvement of no real practical value could be misconstrued as a large improvement.
As an example, if the benchmark has an KGE score of 0.999 and its counterpart 0.9995, then a skill score of 0.5 will
be reported.

# References
1. Knoben, W.J.M., Freer, J.E., Woods, R.A., 2019. 
    Technical note: Inherent benchmark or not? Comparing Nash-Sutcliffe and Kling-Gupta efficiency scores (preprint). 
    Catchment hydrology/Modelling approaches. 
    https://doi.org/10.5194/hess-2019-327

2. Towner, J., Cloke, H.L., Zsoter, E., Flamig, Z., Hoch, J.M., Bazo, J., Coughlan de Perez, E., Stephens, E.M., 2019. 
    Assessing the performance of global hydrological models for capturing peak river flows in the Amazon basin. 
    Hydrology and Earth System Sciences 23, 3057–3080. 
    https://doi.org/10.5194/hess-23-3057-2019
"""
function skill_score(model_score, benchmark_score)
    return (model_score - benchmark_score) / (1.0 - benchmark_score)
end
