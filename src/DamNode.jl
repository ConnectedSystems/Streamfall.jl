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


Base.@kwdef mutable struct DamNode{A <: Union{Param, Real}} <: NetworkNode
    @network_node

    max_storage::Float64
    storage_coef::A = Param(0.5, bounds=(0.00001, 10.0))


    # Function to use to calculate dam level from storage volume
    calc_dam_level::Function = c_dam_level

    # Function to use to calculate dam surface area from storage volume
    calc_dam_area::Function = c_dam_area

    # Function to use to calculate dam discharge from storage volume
    calc_dam_discharge::Function = c_dam_discharge

    # Function to calculate outflow from dam
    calc_dam_outflow::Function = c_dam_outflow

    storage::Array{Float64} = [0.0]

    effective_rainfall::Array{Float64} = []
    et::Array{Float64} = []
    inflow::Array{Float64} = []
    dam_area::Array{Float64} = []

    level::Array{Float64} = []
    discharge::Array{Float64} = []
    outflow::Array{Float64} = []

end


function DamNode(
    name::String,
    area::Float64,
    max_storage::Float64,
    storage_coef::Float64,
    initial_storage::Float64,
    calc_dam_level::Function,
    calc_dam_area::Function,
    calc_dam_discharge::Function,
    calc_dam_outflow::Function)
    return DamNode(name, area, max_storage, storage_coef,
                   calc_dam_level, calc_dam_area, calc_dam_discharge, calc_dam_outflow,
                   [initial_storage], [], [], [], [], [], [], [])
end


"""
    DamNode(node_id::String, spec::Dict)

Create DamNode from a given specification.
"""
function DamNode(node_id::String, spec::Dict)
    n = DamNode{Param}(; node_id=node_id, area=spec["area"], route=false,
                       max_storage=spec["max_storage"])

    node_params = spec["parameters"]
    for (k, p) in node_params
        s = Symbol(k)
        if p isa String
            p = eval(Meta.parse(p))
        end

        try
            if k == "initial_storage"
                setfield!(n, :storage, [p])
            else
                f = getfield(n, s)
                setfield!(n, s, Param(p, bounds=f.bounds))
            end
        catch err
            setfield!(n, s, p)
        end
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


function update_state(node::DamNode, storage, rainfall, et, area, discharge, outflow)
    push!(node.storage, storage)
    push!(node.effective_rainfall, rainfall)
    push!(node.et, et)
    push!(node.level, node.calc_dam_level(storage))
    push!(node.dam_area, area)
    push!(node.discharge, discharge)
    push!(node.outflow, outflow)

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


"""
Calculate outflow for the dam node for a single time step.

# Parameters
- node : DamNode
- rain : rainfall in mm
- et : evapotranspiration data in mm
- irrig_ext : irrigation extractions 
- extractions : extraction data in ML
- gw_flux : groundwater interaction

# Returns
- outflow from dam
"""
function run_node!(node::DamNode, 
                   rain::Float64,
                   et::Float64,
                   inflow::Float64,
                   extractions::Float64,
                   gw_flux::Float64=0.0)

    volume = storage(node)
    dam_area = node.calc_dam_area(volume)
    discharge = node.calc_dam_discharge(volume, node.max_storage)

    updated_store = update_volume(volume, inflow, gw_flux, rain, et,
                                  dam_area, extractions, discharge, node.max_storage)
    outflow = node.calc_dam_outflow(discharge, extractions)

    update_state(node, updated_store, rain, et, dam_area, discharge, outflow)

    return outflow, level(node)
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
    param_info(node::DamNode; kwargs...)::Tuple

Extract node parameter values and bounds.
"""
function param_info(node::DamNode; kwargs...)::Tuple
    tmp = Model(node)
    values = collect(tmp.val)
    bounds = collect(tmp.bounds)
    
    return values, bounds
end


"""
    update_params!(node::DamNode, storage_coef::Float64)::Nothing
"""
function update_params!(node::DamNode, storage_coef::Float64)::Nothing
    node.storage_coef = Param(storage_coef, bounds=node.storage_coef.bounds)
    return nothing
end
