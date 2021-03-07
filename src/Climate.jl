using DataFrames


mutable struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    # t_id::String
end

function subcatchment_data(node::NetworkNode, climate::Climate)
    node_id = node.node_id
    cols = filter(x -> occursin(node_id, string(x)), propertynames(climate))

    return cols
end

function climate_values(node::NetworkNode, climate::Climate, timestep::Union{Nothing, Int}=nothing)
    node_id = node.node_id

    data = climate.climate_data

    rain_col = filter(x -> occursin(node_id, string(x))
                        & occursin(climate.rainfall_id, string(x)),
                        propertynames(data))
    et_col = filter(x -> occursin(node_id, string(x))
                      & occursin(climate.et_id, string(x)),
                      propertynames(data))

    if isempty(rain_col) | isempty(et_col)
        return (missing, missing)
    end

    if isnothing(timestep)
        return select(data, vcat(rain_col, et_col))
    end

    return select(data, vcat(rain_col, et_col))[timestep, :]
end

function sim_length(climate::Climate)
    return nrow(climate.climate_data)
end

