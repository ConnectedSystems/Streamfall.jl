using DataFrames


struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    # t_id::Union{String, Nothing} = nothing
end


"""
    subcatchment_data(node::NetworkNode, climate::Climate)

Extract data for a given node from climate object.
"""
function subcatchment_data(node::NetworkNode, climate::Climate)::DataFrame
    node_id = node.node_id
    data = climate.climate_data
    cols = filter(x -> occursin(node_id, string(x)), names(data))

    return data[:, cols]
end


"""
    climate_values(node::NetworkNode, climate::Climate, timestep::Union{Nothing, Int}=nothing)

Extract climate related data for a given time step.
"""
function climate_values(node::NetworkNode, climate::Climate, 
                        timestep::Union{Nothing, Int}=nothing)
    node_id = node.node_id

    data = climate.climate_data

    rain_col = filter(x -> occursin(node_id, x)
                        & occursin(climate.rainfall_id, x),
                        names(data))[1]
    et_col = filter(x -> occursin(node_id, x)
                      & occursin(climate.et_id, x),
                      names(data))[1]

    if isempty(rain_col) | isempty(et_col)
        return (missing, missing)
    end

    if isnothing(timestep)
        return select(data, [rain_col, et_col])
    end

    return data[timestep, [rain_col, et_col]]
end


"""
    sim_length(climate::Climate)::Int64

Simulation length is dependent on available climate data.
"""
function sim_length(climate::Climate)::Int64
    return nrow(climate.climate_data)
end

