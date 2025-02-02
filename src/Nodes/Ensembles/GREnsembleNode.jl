Base.@kwdef mutable struct GREnsembleNode{N<:NetworkNode, P, A<:Real} <: EnsembleNode
    name::String
    area::A

    instances::Array{N} = NetworkNode[]

    # GRC method
    comb_method::Function = grc_combine

    outflow::Array{A} = []

    obj_func::Function = obj_func
end

function GREnsembleNode(nodes::Vector{<:NetworkNode})::GREnsembleNode
    n1 = nodes[1]
    gr_node = GREnsembleNode{NetworkNode, Param, Float64}(;
        name=n1.name,
        area=n1.area,
        instances=nodes
    )

    return gr_node
end

function create_node(
    node::Type{<:GREnsembleNode},
    nodes::Vector{<:NetworkNode};
    kwargs...
)
    return GREnsembleNode(nodes; kwargs...)
end

function grc_weights(X::Matrix{T}, y::Vector{T}) where T<:Real
    # Add constant term for bias correction
    X_aug = hcat(ones(size(X,1)), X)

    # Solve normal equations: β = (X'X)^(-1)X'y
    β = inv(X_aug' * X_aug) * X_aug' * y

    # Split bias term and weights
    bias = β[1]
    weights = β[2:end]

    return weights, bias
end

function grc_combine(X::Matrix{T}, weights::Vector{T}, bias::T) where T<:Real
    # Apply weights and bias correction
    return X * weights .+ bias
end

function calibrate_instances!(
    ensemble::GREnsembleNode,
    climate::Climate,
    calib_data::DataFrame,
    metric::Union{F,AbstractDict{String,F}};
    kwargs...
) where {F}

    # Calibrate individual instances first
    for node in ensemble.instances
        calibrate!(node, climate, calib_data, metric; kwargs...)
    end

    # Then determine GRC
    return calibrate!(ensemble, climate, calib_data, metric; kwargs...)
end

function calibrate!(
    ensemble::GREnsembleNode,
    climate::Climate,
    calib_data::DataFrame,
    metric::Union{C,AbstractDict{String,C}};  # Unused, added to maintain consistent interface
    kwargs...
) where {C<:Function}

    for inst in ensemble.instances
        run_node!(inst, climate)
    end

    X = Matrix(hcat([m.outflow for m in ensemble.instances]...))
    weights, bias = grc_weights(X, calib_data[:, ensemble.name])

    ensemble.comb_method = (X) -> grc_combine(hcat(X...), weights, bias)

    return nothing
end

function run_node!(ensemble::GREnsembleNode, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    for inst in ensemble.instances
        run_node!(inst, climate; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    X = hcat([m.outflow for m in ensemble.instances]...)'

    ensemble.outflow = ensemble.comb_method([inst.outflow for inst in ensemble.instances])
end

function run_timestep!(node::WeightedEnsembleNode, rain, et, ts; inflow=0.0, extraction=0.0, exchange=0.0)
    for inst in node.instances
        run_timestep!(inst, rain, et, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    node.outflow[ts] = node.comb_method([inst.outflow[ts] for inst in node.instances])
end
