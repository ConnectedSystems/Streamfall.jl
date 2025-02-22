using ModelParameters


function c_dam_level(volume)
    return 156.8 + 0.9463 * volume^0.2922
end


function c_dam_area(volume)
    return 0.0021 * volume^0.762
end


function c_dam_discharge(volume, max_storage)
    discharge = 0.0
    if volume > max_storage
        discharge = 0.001492 * (volume - max_storage)^1.5280
    end

    return max(0.0, discharge)
end


function c_dam_outflow(discharge, irrigation_extraction)
    return discharge + irrigation_extraction
end


Base.@kwdef mutable struct DamNode{P,A<:AbstractFloat} <: NetworkNode
    name::String
    area::A

    max_storage::A
    storage_coef::P = Param(0.5, bounds=(0.00001, 10.0))

    # Dam storage volume to level
    calc_dam_level::Function = c_dam_level

    # Dam volume to surface area
    calc_dam_area::Function = c_dam_area

    # Dam storage volume to dam discharge
    calc_dam_discharge::Function = c_dam_discharge

    # Dam outflow
    calc_dam_outflow::Function = c_dam_outflow

    storage::Array{A} = [0.0]

    effective_rainfall::Array{A} = []
    et::Array{A} = []
    inflow::Array{A} = []
    dam_area::Array{A} = []

    level::Array{A} = []
    discharge::Array{A} = []
    outflow::Array{A} = []

    obj_func::Function = dependent_obj_func

    function_repr = OrderedDict()
end


function DamNode(
    name::String,
    area::F,
    max_storage::F,
    storage_coef::F,
    initial_storage::F,
    calc_dam_level::Function,
    calc_dam_area::Function,
    calc_dam_discharge::Function,
    calc_dam_outflow::Function
) where {F<:Float64}
    return DamNode(name, area, max_storage, storage_coef,
        calc_dam_level, calc_dam_area, calc_dam_discharge, calc_dam_outflow,
        F[initial_storage], F[], F[], F[], F[], F[], F[], F[])
end

"""
    DamNode(name::String, spec::AbstractDict)

Create DamNode from a given specification.
"""
function DamNode(name::String, spec::AbstractDict)
    n = DamNode{Param,Float64}(;
        name=name,
        area=spec["area"],
        max_storage=spec["max_storage"]
    )
    n.storage = [spec["initial_storage"]]

    node_params = spec["parameters"]
    for (k, p) in node_params
        s = Symbol(k)
        if p isa String
            p = eval(Meta.parse(p))
        end

        f = getfield(n, s)
        setfield!(n, s, Param(p, bounds=f.bounds))
    end

    for (k, p) in spec["functions"]
        s = Symbol(k)
        f = getfield(n, s)

        n.function_repr[k] = p

        p = eval(Meta.parse(p))
        setfield!(n, s, p)
    end

    return n
end


function level(node::DamNode)
    return last(node.level)
end


function area(node::DamNode)
    return last(node.dam_area)
end


function discharge(node::DamNode)
    return last(node.discharge)
end


function storage(node::DamNode)
    return last(node.storage)
end

function prep_state!(node::DamNode, timesteps::Int64)
    resize!(node.storage, timesteps + 1)
    node.storage[2:end] .= 0.0

    node.effective_rainfall = zeros(timesteps)
    node.et = zeros(timesteps)
    node.inflow = zeros(timesteps)
    node.dam_area = zeros(timesteps)

    node.level = zeros(timesteps)
    node.discharge = zeros(timesteps)
    node.outflow = zeros(timesteps)
end


function update_state!(node::DamNode, storage, rainfall, et, area, discharge, outflow)
    push!(node.storage, storage)
    push!(node.effective_rainfall, rainfall)
    push!(node.et, et)
    push!(node.level, node.calc_dam_level(storage))
    push!(node.dam_area, area)
    push!(node.discharge, discharge)
    push!(node.outflow, outflow)

    return nothing
end
function update_state!(node::DamNode, ts::Int64, storage, rainfall, et, area, discharge, outflow)
    node.storage[ts+1] = storage

    node.effective_rainfall[ts] = rainfall
    node.et[ts] = et
    node.level[ts] = node.calc_dam_level(storage)
    node.dam_area[ts] = area
    node.discharge[ts] = discharge
    node.outflow[ts] = outflow

    return nothing
end


"""
    update_volume(volume, node_inflow, gamma, rain, evap, area, extractions, discharge, max_store)::Float64

Update dam volume for timestep.

# Arguments
- volume : current water volume in ML
- node_inflow : inflow from previous node in ML
- gamma : groundwater exchange (positive is gain from gw flow, negative is loss to infiltration)
- rain : rainfall input
- evap : evaporation loss
- infiltration : infiltration loss
- area : dam surface area in square kilometers
- extractions : water extraction from dam in ML
- discharge : discharge from dam in ML
- max_store : maximum dam storage in ML

# Returns
volume of water stored in dam
"""
function update_volume(volume, node_inflow, gamma, rain, evap, area, extractions, discharge, max_store)::Float64

    vol = volume + (node_inflow + gamma) + (rain - evap) * area - extractions - discharge

    return max(0.0, min(max_store, vol))
end


function run_node!(node::DamNode, climate::Climate;
    inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = sim_length(climate)
    prep_state!(node, timesteps)

    for ts in 1:timesteps
        run_node!(node, climate, ts;
            inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return nothing
end


"""
    run_node!(
        node::DamNode, climate::Climate, ts::Int;
        inflow=nothing, extraction=nothing, exchange=nothing
    )::Nothing

Run a specific node for a specified time step.

# Arguments
- `node` : DamNode
- `climate` : Climate dataset
- `ts` : Current time step
- `inflow` : Time series of inflows from any upstream node.
- `extraction` : Time series of water orders (expects column of `_releases`)
- `exchange` : Time series of groundwater flux
"""
function run_node!(
    node::DamNode, climate::Climate, ts::Int;
    inflow=nothing, extraction=nothing, exchange=nothing
)::Nothing
    node_name = node.name
    rain, et = climate_values(node, climate, ts)
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)
    current_vol = node.storage[ts]

    run_node!(node, ts, rain, et, current_vol, in_flow, wo, ex)

    return nothing
end


"""
    run_node!(
        node::DamNode,
        ts::Int64,
        rain::Float64,
        et::Float64,
        volume::Float64,
        inflow::Float64,
        extractions::Float64,
        gw_flux::Float64
    )

Calculate outflow for the dam node for a single time step.

# Arguments
- `node` : DamNode
- `rain` : rainfall in mm
- `et` : evapotranspiration data in mm
- `irrig_ext` : irrigation extractions
- `extractions` : extraction data in ML
- `gw_flux` : groundwater interaction

# Returns
Outflow from dam
"""
function run_node!(
    node::DamNode,
    ts::Int64,
    rain::Float64,
    et::Float64,
    volume::Float64,
    inflow::Float64,
    extractions::Float64,
    gw_flux::Float64
)
    dam_area = node.calc_dam_area(volume)
    discharge = node.calc_dam_discharge(volume, node.max_storage)

    updated_store = update_volume(volume, inflow, gw_flux, rain, et,
        dam_area, extractions, discharge, node.max_storage)
    outflow = node.calc_dam_outflow(discharge, extractions)

    update_state!(node, ts, updated_store, rain, et, dam_area, discharge, outflow)

    return nothing
end


function reset!(node::DamNode)::Nothing
    node.storage = [node.storage[1]]
    node.effective_rainfall = []
    node.et = []
    node.level = []
    node.dam_area = []

    node.discharge = []
    node.outflow = []

    return nothing
end


"""
    update_params!(node::DamNode, storage_coef::Float64)::Nothing

Method to update `DamNode` specific parameters.

# Arguments
- `node` : DamNode
- `storage_coef` : Storage coefficient value
"""
function update_params!(node::DamNode, storage_coef::Float64)::Nothing
    node.storage_coef = Param(storage_coef, bounds=node.storage_coef.bounds)
    return nothing
end

"""
    extract_spec!(node::DamNode, spec::AbstractDict)::Nothing

Extract dam-specific values.
"""
function extract_spec!(node::DamNode, spec::AbstractDict)::Nothing
    spec["initial_storage"] = node.storage[1]
    spec["max_storage"] = node.max_storage
    spec["functions"] = OrderedDict()

    for (f, r) in node.function_repr
        spec["functions"][f] = r
    end

    return nothing
end
