using CSV
using DataFrames


struct Climate
    climate_data::DataFrame
    rainfall_id::String
    et_id::String
    t_id::String
end
function Climate(data::DataFrame, p_id, et_id)
    return Climate(data, p_id, et_id, "_T")
end

function Climate(file_path::String, p_id, et_id; t_id="_T")
    climate_data = CSV.read(file_path, DataFrame; comment="#")

    return Climate(climate_data, p_id, et_id, t_id)
end

"""
    extract_flow(
        data::DataFrame, gauge_id::String; suffix::String="_Q"
    )::DataFrame

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
    data::DataFrame, gauge_id::String; suffix::String="_Q"
)::DataFrame
    target = data[:, ["Date", gauge_id * suffix]]
    try
        target[!, gauge_id*suffix] = convert.(Float64, target[!, gauge_id*suffix])
    catch
        target[!, gauge_id*suffix] = convert.(Union{Float64,Missing}, target[!, gauge_id*suffix])
    end

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
- `T_suffix` : Suffix used to indicate Temperature (default: "_PET")

# Returns
Climate
"""
@inline function extract_climate(
    data::DataFrame; P_suffix::String="_P", PET_suffix::String="_PET", T_suffix::String="_T"
)::Climate
    return Climate(data, P_suffix, PET_suffix, T_suffix)
end

"""
    rainfall_data(node::NetworkNode, climate::Climate)::DataFrame

Extract rainfall data for a given node.
"""
function rainfall_data(node::NetworkNode, climate::Climate)::DataFrame
    data = climate.climate_data
    rain_col = filter(x -> occursin(node.name, x)
                           &
                           occursin(climate.rainfall_id, x),
        names(data))[1]

    return data[:, rain_col]
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
    climate_values(node::NetworkNode, climate::Climate, timestep::Int)

Extract climate related data for a given time step.
"""
function climate_values(node::NetworkNode, climate::Climate, timestep::Int)
    node_name::String = node.name
    data::DataFrame = climate.climate_data

    # TODO : Catch instances where data is not found (raises BoundsError)
    rain_col = filter(x -> occursin(node_name, x)
                           &
                           occursin(climate.rainfall_id, x),
        names(data))[1]
    et_col = filter(x -> occursin(node_name, x)
                         &
                         occursin(climate.et_id, x),
        names(data))[1]
    t_col = try
        filter(x -> occursin(node_name, x)
                    &
                    occursin(climate.t_id, x),
            names(data))[1]
    catch err
        if !(err isa BoundsError)
            rethrow(err)
        end

        []
    end

    if isempty(rain_col) | isempty(et_col)
        throw(ArgumentError("No climate data found for $(node_name) at time step: $(timestep)"))
    end

    if isempty(t_col)
        sel = [rain_col, et_col]
    else
        sel = [rain_col, et_col, t_col]
    end

    return data[timestep, sel]
end


"""
    climate_values(node::NetworkNode, climate::Climate)

Extract climate related data.
"""
function climate_values(node::NetworkNode, climate::Climate)
    data::DataFrame = climate.climate_data

    # TODO : Catch instances where data is not found (raises BoundsError)
    rain_col = filter(x -> occursin(node.name, x)
                           &
                           occursin(climate.rainfall_id, x),
        names(data))[1]
    et_col = filter(x -> occursin(node.name, x)
                         &
                         occursin(climate.et_id, x),
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

Base.show(io::IO, ::MIME"text/plain", c::Climate) = show(io, c)
function Base.show(io::IO, c::Climate)
    ntype_name = nameof(typeof(c))

    col_names = names(c.climate_data)

    dt_cols = ["year", "month", "day", "Date"]
    ignore_cols = dt_cols[[in(t, col_names) for t in dt_cols]]
    tgt_col_names = names(c.climate_data[:, Not(ignore_cols...)])

    P_cols = occursin.(c.rainfall_id, col_names)
    ET_cols = occursin.(c.et_id, col_names)
    T_cols = occursin.("_T", col_names)

    ts_diff = diff(c.climate_data.Date)
    is_contiguous = all(ts_diff .== ts_diff[1])  # False if there are any breaks in the data

    first_period = first(c.climate_data[:, "Date"])
    last_period = last(c.climate_data[:, "Date"])
    approx_years = round((last_period - first_period).value / 365.25; digits=2)
    println(io, "Climate dataset\n")
    println(io, "Number of non-date columns: $(length(tgt_col_names))")
    println(io, "    P identifier: \"$(c.rainfall_id)\"")
    println(io, "    ET identifier: \"$(c.et_id)\"")
    println(io, "    T identifier: \"$(c.t_id)\"")
    print(io, "\n")
    println(io, "Represented time period: $(first_period) - $(last_period)")
    println(io, "    Number of observations: $(length(c))")
    println(io, "    Approx. number of years: $(approx_years)")
    println(io, "    Is contiguous: $(is_contiguous)")
    print(io, "\n")
    println(io, "Preview:")
    println(io, first(c.climate_data, 10))

    print(io, "\n")
end