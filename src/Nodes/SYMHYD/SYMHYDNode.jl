using Parameters
using ModelParameters


const SYMHYD_SOIL_ET_CONST = 10.0


"""
"""
Base.@kwdef mutable struct SYMHYDNode{P, A<:AbstractFloat} <: NetworkNode
    name::String
    area::A

    # parameters
    baseflow_coef::P = Param(0.5, bounds=(0.0, 1.0))
    impervious_threshold::P = Param(2.5, bounds=(0.0, 5.0))  # mm
    infiltration_coef::P = Param(200.0, bounds=(0.0, 400.0))
    infiltration_shape::P = Param(5.0, bounds=(0.0, 10.0))
    interflow_coef::P = Param(0.5, bounds=(0.0, 1.0))
    pervious_fraction::P = Param(0.5, bounds=(0.0, 1.0))
    risc::P = Param(2.5, bounds=(0.0, 5.0))  # rainfall interception store capacity (mm)
    recharge_coef::P = Param(0.5, bounds=(0.0, 1.0))
    smsc::P = Param(250.0, bounds=(1.0, 500.0))  # Soil Moisture Store Capacity (mm)

    # stores
    sm_store::Array{A} = [0.0]
    gw_store::Array{A} = [0.0]
    total_store::Array{A} = [0.0]

    # outputs
    outflow::Array{A} = []  # mm
    baseflow::Array{A} = []  # mm
    quickflow::Array{A} = []  # mm

    obj_func::Function = obj_func
end


"""

Create node from spec.
"""
function SYMHYDNode(name::String, spec::AbstractDict)
    n = create_node(SYMHYDNode, name, spec["area"])
    node_params = spec["parameters"]
    node_params["sm_store"] = [node_params["initial_sm_store"]]
    node_params["gw_store"] = [node_params["initial_gw_store"]]
    node_params["total_store"] = [node_params["initial_total_store"]]

    remove_from_spec = (x) -> delete!(node_params, x)
    map(remove_from_spec, ["initial_sm_store", "initial_gw_store", "initial_total_store"])

    for (k, p) in node_params
        s = Symbol(k)
        if p isa String
            p = eval(Meta.parse(p))
        end

        try
            f = getfield(n, s)
            setfield!(n, s, Param(p, bounds=f.bounds))
        catch err
            msg = sprint(showerror, err, catch_backtrace())
            if occursin("no field bounds", string(msg))
                setfield!(n, s, p)
            else
                throw(err)
            end
        end
    end

    return n
end

function prep_state!(node::SYMHYDNode, sim_length::Int64)::Nothing
    resize!(node.sm_store, sim_length+1)
    resize!(node.gw_store, sim_length+1)
    resize!(node.total_store, sim_length+1)
    node.sm_store[2:end] .= 0.0
    node.gw_store[2:end] .= 0.0
    node.total_store[2:end] .= 0.0

    node.outflow = fill(0.0, sim_length)
    node.baseflow = fill(0.0, sim_length)
    node.quickflow = fill(0.0, sim_length)

    return nothing
end


"""

Run SYMHYD for a given time step
"""
function run_timestep!(node::SYMHYDNode, climate::Climate, ts::Int; inflow=nothing, extraction=extraction, exchange=nothing)::AbstractFloat
    P, E = climate_values(node, climate, ts)

    return run_timestep!(node, P, E, ts; inflow=inflow, extraction=extraction, exchange=exchange)
end

function run_timestep!(
    node::SYMHYDNode,
    rain::F,
    et::F,
    ts::Int;
    inflow=nothing,
    extraction=nothing,
    exchange=nothing
)::F where {F<:AbstractFloat}
    sm_store, gw_store, total_store, total_runoff, baseflow, event_runoff = run_symhyd(node, rain, et, ts)

    node_name = node.name
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)

    if !isnothing(inflow)
        total_runoff = total_runoff + in_flow + ex - wo
    end

    area::F = node.area
    update_state!(node, ts, sm_store, gw_store, total_store, total_runoff * area, baseflow * area, event_runoff * area)

    return total_runoff
end


function update_state!(node::SYMHYDNode, sm_store, gw_store, total_store, outflow, baseflow, quickflow)
    push!(node.sm_store, sm_store)
    push!(node.gw_store, gw_store)
    push!(node.total_store, total_store)
    push!(node.outflow, outflow)
    push!(node.baseflow, baseflow)
    push!(node.quickflow, quickflow)
end
function update_state!(node::SYMHYDNode, ts::Int64, sm_store, gw_store, total_store, outflow, baseflow, quickflow)
    node.sm_store[ts+1] = sm_store
    node.gw_store[ts+1] = gw_store
    node.total_store[ts+1] = total_store
    node.outflow[ts] = outflow
    node.baseflow[ts] = baseflow
    node.quickflow[ts] = quickflow
end


function update_params!(node::SYMHYDNode, baseflow_coef::Float64, impervious_threshold::Float64,
                        infiltration_coef::Float64,
                        infiltration_shape::Float64,
                        interflow_coef::Float64,
                        pervious_fraction::Float64,
                        risc::Float64,
                        recharge_coef::Float64,
                        smsc::Float64)::Nothing
    node.baseflow_coef = Param(baseflow_coef, bounds=node.baseflow_coef.bounds)
    node.impervious_threshold = Param(impervious_threshold, bounds=node.impervious_threshold.bounds)
    node.infiltration_coef = Param(infiltration_coef, bounds=node.infiltration_coef.bounds)
    node.infiltration_shape = Param(infiltration_shape, bounds=node.infiltration_shape.bounds)
    node.interflow_coef = Param(interflow_coef, bounds=node.interflow_coef.bounds)
    node.pervious_fraction = Param(pervious_fraction, bounds=node.pervious_fraction.bounds)
    node.risc = Param(risc, bounds=node.risc.bounds)
    node.recharge_coef = Param(recharge_coef, bounds=node.recharge_coef.bounds)
    node.smsc = Param(smsc, bounds=node.smsc.bounds)

    return nothing
end


"""
    reset!(node::SYMHYDNode)::Nothing

Reset node. Clears all states back to their initial values.
"""
function reset!(node::SYMHYDNode)::Nothing
    # stores
    node.sm_store = [node.sm_store[1]]
    node.gw_store = [node.gw_store[1]]
    node.total_store = [node.total_store[1]]

    # outputs
    node.outflow = []
    node.baseflow = []
    node.quickflow = []

    return nothing
end


"""

Run SYMHYD for a single time step with given inputs and state variables.
"""
function run_symhyd(node::SYMHYDNode, P::F, ET::F, ts::Int64)::NTuple{6,F} where {F<:Float64}

    sm_store::F = node.sm_store[ts]
    gw_store::F = node.gw_store[ts]
    total_store::F = node.total_store[ts]

    pervious_incident::F = P
	impervious_incident::F = P

    impervious_ET::F = min(node.impervious_threshold, impervious_incident)
    impervious_runoff::F = impervious_incident - impervious_ET

    interception_ET::F = min(pervious_incident, min(ET, node.risc))
    throughfall::F = pervious_incident - impervious_ET

    smsc::F = node.smsc.val
    sm_fraction::F = sm_store / smsc

    infiltration_capacity::F = node.infiltration_coef.val::F * exp(-node.infiltration_shape.val::F * sm_fraction)::F
    infiltration::F = min(throughfall, infiltration_capacity)
    infiltration_Xs_runoff::F = throughfall - infiltration

    interflow_runoff::F = node.interflow_coef.val::F * sm_fraction * infiltration
    infiltration_after_interflow::F = infiltration - interflow_runoff

    recharge::F = node.recharge_coef.val::F * sm_fraction * infiltration_after_interflow
    soil_input::F = infiltration_after_interflow - recharge

    # Update states...
    sm_store += soil_input
    sm_fraction = sm_store / smsc
    gw_store += recharge

    if sm_fraction > 1.0
        gw_store += sm_store - smsc
        sm_store = smsc
        sm_fraction = 1.0
    end

    baseflow_runoff::F = node.baseflow_coef.val::F * gw_store
    gw_store -= baseflow_runoff

    soil_ET::F = min(sm_store, min(ET - interception_ET, sm_fraction*SYMHYD_SOIL_ET_CONST))
    sm_store -= soil_ET

    pervious_frac::F = node.pervious_fraction.val::F
    total_store = (sm_store + gw_store) * pervious_frac

    event_runoff::F = (1.0 - pervious_frac) * impervious_runoff + pervious_frac * (infiltration_Xs_runoff + interflow_runoff)
    total_runoff::F = event_runoff + pervious_frac * baseflow_runoff
    baseflow::F = baseflow_runoff * pervious_frac

	# values for time step...
	return sm_store, gw_store, total_store, total_runoff, baseflow, event_runoff
end
