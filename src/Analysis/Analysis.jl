module Analysis
include("uncertainty.jl")
end


"""
    gmd(vals::AbstractVector{<:Real})::Float64
    gmd(vals::AbstractMatrix{<:Real})

Gini's Mean Difference.

The absolute mean of all pairwise distances between elements in a given set.

# References
1. La Haye, R., & Zizler, P. (2019). 
   The Gini mean difference and variance. 
   METRON, 77(1), 43-52. 
   https://doi.org/10.1007/s40300-019-00149-2

2. Yitzhaki, S. (2003). 
   Gini's Mean difference: A superior measure of variability for non-normal
     distributions. 
   Metron - International Journal of Statistics, LXI(2), 285-316.
   https://ideas.repec.org/a/mtn/ancoec/030208.html

3. Kashif, M., Aslam, M., Al-Marshadi, A. H., & Jun, C.-H. (2016).
   Capability Indices for Non-Normal Distribution Using Gini's Mean Difference as Measure of Variability. 
   IEEE Access, 4, 7322-7330.
   https://doi.org/10.1109/ACCESS.2016.2620241
"""
function gmd(vals::AbstractVector{<:Real})::Float64
    n = length(vals)
    sv = sort(vals)
    return (2 / (n * (n - 1))) .* sum(([((2 * i) - n - 1) * sv[i] for i in 1:n]))
end
function gmd(vals::AbstractMatrix{<:Real})
    return gmd.(eachcol(vals))
end
