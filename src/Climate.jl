using DataFrames


struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    # t_id::Union{String, Nothing} = nothing
end

"""
Extract streamflow data from file.

Streamflow (Q) column is identified the Gauge ID.

Flow data is identified with the suffix `_Q` by default
e.g., ("000001_Q")

# Arguments
- `data` : Observation data
- `gauge_id` : Gauge/Node ID
- `suffix` : Suffix used to indicate flow data (default: "_Q")

# Returns
DataFrame of observations for selected gauge.
"""
@inline function extract_flow(
    data::DataFrame, gauge_id::String, suffix::String="_Q"
)::DataFrame
    target = data[:, ["Date", gauge_id * suffix]]
    rename!(target, gauge_id * suffix => gauge_id)

    return target
end

"""
Create a climate dataset of Precipitation (P) and Potential Evapotranspiration (PET).
Data for multiple gauges may be defined in a single dataset.

P and PET columns are identified by `_P` and `_PET` suffixes by default.

# Arguments
- `data` : Observation data
- `P_suffix` : Suffix used to indicate precipitation (default: "_P")
- `PET_suffix` : Suffix used to indicate Potential Evapotranspiration (default: "_PET")

# Returns
Climate
"""
@inline function extract_climate(
    data::DataFrame; P_suffix::String="_P", PET_suffix::String="_PET"
)::Climate
    return Climate(data, P_suffix, PET_suffix)
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

Extract climate related data.
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
