using DataFrames


struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    # t_id::Union{String, Nothing} = nothing
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

    # if !isnothing(climate.t_id)
    #     t_col = filter(x -> occursin(node_id, string(x))
    #                   & occursin(climate.t_id, string(x)),
    #                   propertynames(data))
    # end

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

