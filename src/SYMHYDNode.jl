using Parameters
using ModelParameters


const SYMHYD_SOIL_ET_CONST = 10.0


"""
"""
Base.@kwdef mutable struct SYMHYDNode{P} <: NetworkNode
    @network_node

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
    sm_store::Array{Float64} = [0.0]
    gw_store::Array{Float64} = [0.0]
    total_store::Array{Float64} = [0.0]

    # outputs
    outflow::Array{Float64} = []  # mm
    baseflow::Array{Float64} = []  # mm
    quickflow::Array{Float64} = []  # mm
end


"""

Create node from spec.
"""
function SYMHYDNode(name::String, spec::Dict)
    n = SYMHYDNode{Param}(; name=name, area=spec["area"])
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


"""

Run SYMHYD for all time steps within a climate scenario.
"""
function run_node!(node::SYMHYDNode, climate::Climate; inflow=nothing, extraction=nothing, exchange=nothing)
    timesteps = length(climate)
    for ts in 1:timesteps
        run_node!(node, climate, ts; inflow=inflow, extraction=extraction, exchange=exchange)
    end

    return node.outflow
end


"""

Run SYMHYD for a given time step
"""
function run_node!(node::SYMHYDNode, climate::Climate, timestep::Int; inflow=nothing, extraction=extraction, exchange=nothing)
    P, E = climate_values(node, climate, timestep)

    res = run_symhyd(node, P, E)
    sm_store, gw_store, total_store, total_runoff, baseflow, event_runoff = res

    ts = timestep
    node_name = node.name
    wo = timestep_value(ts, node_name, "releases", extraction)
    ex = timestep_value(ts, node_name, "exchange", exchange)
    in_flow = timestep_value(ts, node_name, "inflow", inflow)

    if !isnothing(inflow)
        total_runoff = total_runoff + in_flow + ex - wo
    end

    area = node.area
    update_state!(node, sm_store, gw_store, total_store, total_runoff * area, baseflow * area, event_runoff * area)
end


function update_state!(node::SYMHYDNode, sm_store, gw_store, total_store, outflow, baseflow, quickflow)
    push!(node.sm_store, sm_store)
    push!(node.gw_store, gw_store)
    push!(node.total_store, total_store)
    push!(node.outflow, outflow)
    push!(node.baseflow, baseflow)
    push!(node.quickflow, quickflow)
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
function run_symhyd(node::SYMHYDNode, P, ET)

    sm_store = node.sm_store[end]
    gw_store = node.gw_store[end]
    total_store = node.total_store[end]

    pervious_incident = P
	impervious_incident = P

    impervious_ET = min(node.impervious_threshold, impervious_incident)
    impervious_runoff = impervious_incident - impervious_ET

    interception_ET = min(pervious_incident, min(ET, node.risc))
    throughfall = pervious_incident - impervious_ET

    smsc = node.smsc
    sm_fraction = sm_store / smsc

    infiltration_capacity = node.infiltration_coef * exp(-node.infiltration_shape*sm_fraction)
    infiltration = min(throughfall, infiltration_capacity)
    infiltration_Xs_runoff = throughfall - infiltration

    interflow_runoff = node.interflow_coef * sm_fraction * infiltration
    infiltration_after_interflow = infiltration - interflow_runoff

    recharge = node.recharge_coef * sm_fraction * infiltration_after_interflow
    soil_input = infiltration_after_interflow - recharge

    # Update states...
    sm_store += soil_input
    sm_fraction = sm_store / smsc
    gw_store += recharge

    if sm_fraction > 1.0
        gw_store += sm_store - smsc
        sm_store = smsc
        sm_fraction = 1.0
    end

    baseflow_runoff = node.baseflow_coef * gw_store
    gw_store -= baseflow_runoff

    soil_ET = min(sm_store, min(ET - interception_ET, sm_fraction*SYMHYD_SOIL_ET_CONST))
    sm_store -= soil_ET

    pervious_frac = node.pervious_fraction
    total_store = (sm_store + gw_store) * pervious_frac

    event_runoff = (1.0 - pervious_frac) * impervious_runoff + pervious_frac * (infiltration_Xs_runoff + interflow_runoff)
    total_runoff = event_runoff + pervious_frac * baseflow_runoff
    baseflow = baseflow_runoff * pervious_frac

	# values for time step...
	return sm_store, gw_store, total_store, total_runoff, baseflow, event_runoff
end
