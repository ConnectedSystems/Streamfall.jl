using DataFrames


struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    # t_id::Union{String, Nothing} = nothing
end


"""
    subcatchment_data(node::NetworkNode, climate::Climate)

Extract all data for a given node from climate object.
"""
function subcatchment_data(node::NetworkNode, climate::Climate)::DataFrame
    data = climate.climate_data
    cols = filter(x -> occursin(node.name, string(x)), names(data))

    return data[:, vcat(["Date"], cols)]
end


"""
    rainfall_data(node::NetworkNode, climate::Climate)::DataFrame

Extract rainfall data for a given node.
"""
function rainfall_data(node::NetworkNode, climate::Climate)::DataFrame
    data = climate.climate_data
    rain_col = filter(x -> occursin(node.name, x) 
                            & occursin(climate.rainfall_id, x), 
                            names(data))[1]

    return data[:, rain_col]
end


"""
    climate_values(node::NetworkNode, climate::Climate, timestep::Int)

Extract climate related data for a given time step.
"""
function climate_values(node::NetworkNode, climate::Climate, 
                        timestep::Int)
    node_name::String = node.name
    data::DataFrame = climate.climate_data

    # TODO : Catch instances where data is not found (raises BoundsError)
    rain_col = filter(x -> occursin(node_name, x) 
                        & occursin(climate.rainfall_id, x), 
                        names(data))[1]
    et_col = filter(x -> occursin(node_name, x)
                        & occursin(climate.et_id, x),
                        names(data))[1]

    if isempty(rain_col) | isempty(et_col)
        throw(ArgumentError("No climate data found for $(node_name) at time step: $(timestep)"))
    end

    return data[timestep, [rain_col, et_col]]
end


"""
    climate_values(node::NetworkNode, climate::Climate)

Extract climate related data for a given time step.
"""
function climate_values(node::NetworkNode, climate::Climate)
    data::DataFrame = climate.climate_data

    # TODO : Catch instances where data is not found (raises BoundsError)
    rain_col = filter(x -> occursin(node.name, x) 
                            & occursin(climate.rainfall_id, x), 
                            names(data))[1]
    et_col = filter(x -> occursin(node.name, x)
                            & occursin(climate.et_id, x),
                            names(data))[1]

    if isempty(rain_col) | isempty(et_col)
        throw(ArgumentError("No climate data found for $(node.name) at time step: $(timestep)"))
    end

    return select(data, [rain_col, et_col])
end


"""
    sim_length(climate::Climate)::Int64

Simulation length is dependent on available climate data.
"""
function sim_length(climate::Climate)::Int64
    return nrow(climate.climate_data)
end


function timesteps(climate::Climate)::Array
    return climate.climate_data[:, :Date]
end


Base.length(climate::Climate) = nrow(climate.climate_data)
