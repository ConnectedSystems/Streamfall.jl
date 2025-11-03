using Parameters
using ModelParameters

include("./components/gw_flow.jl")


"""
    IHACRESBilinearNodeGW

IHACRES node with explicit groundwater storage module (IHACRES_GW).

Extends standard IHACRES with spatially-explicit groundwater component that:
- Tracks groundwater storage volume (ML)
- Calculates baseflow from storage-discharge relationship
- Allows groundwater extraction and external fluxes
- Converts storage to water table depth for calibration to bore observations

Based on Ivkovic et al. (2005) IHACRES_GW formulation.

# Fields
## Identifiers
- `name::String` : Node identifier
- `area::Float64` : Catchment area (km²) - used for runoff calculations
- `aquifer_area::Float64` : Aquifer area (km²) - used for GW depth conversions (defaults to catchment area)

## CMD Parameters (same as standard IHACRES)
- `d::P` : CMD threshold (mm), controls runoff generation
- `d2::P` : Scaling factor for second threshold
- `e::P` : ET conversion factor
- `f::P` : Plant stress threshold factor
- `alpha::P` : Effective rainfall scaling factor

## Quick Flow Parameter
- `a::P` : Quickflow coefficient (1/τ_q)

## Groundwater Parameters (NEW)
- `tau_s::P` : Slow flow time constant (days), controls baseflow recession
- `vs::P` : Recharge fraction (0-1), proportion of effective rainfall to GW
- `L::P` : Constant loss term (ML/day), positive=loss, negative=gain
- `specific_yield::P` : Aquifer specific yield (dimensionless, typically 0.1-0.3)

## Elevation Parameters (for depth conversion)
- `gauge_elevation::Float64` : Stream bed elevation at gauge (mAHD)
- `bore_ground_elevation::Float64` : Ground surface elevation at bore (mAHD)

## State Variables
- `storage::Array{Float64}` : CMD (mm), length = timesteps + 1
- `quick_store::Array{Float64}` : Quickflow storage (ML), length = timesteps + 1
- `gw_storage::Array{Float64}` : Groundwater storage (ML), length = timesteps + 1
- `effective_rainfall::Array{Float64}` : Effective rainfall (mm/day), length = timesteps
- `et::Array{Float64}` : Evapotranspiration (mm/day), length = timesteps
- `quickflow::Array{Float64}` : Quickflow discharge (ML/day), length = timesteps
- `baseflow::Array{Float64}` : Baseflow discharge (ML/day), length = timesteps
- `outflow::Array{Float64}` : Total outflow (ML/day), length = timesteps
- `inflow::Array{Float64}` : Inflow from upstream (ML/day), length = timesteps

## Calibration
- `obj_func::Function` : Objective function for calibration
"""
Base.@kwdef mutable struct IHACRESBilinearNodeGW{P,A<:AbstractFloat} <: IHACRESNode
    const name::String
    const area::A
    aquifer_area::A = area  # Defaults to catchment area if not specified

    # CMD parameters (same as standard IHACRES)
    d::P = Param(200.0, bounds=(10.0, 200.0), desc="CMD threshold")
    d2::P = Param(2.0, bounds=(0.0001, 10.0), desc="Scaling factor")
    e::P = Param(1.0, bounds=(0.9, 1.0), desc="ET conversion factor")
    f::P = Param(0.8, bounds=(0.01, 2.5), desc="Plant stress threshold")
    alpha::P = Param(0.1, bounds=(0.1, 1 - 1 / 10^9), desc="Effective rainfall scaling")

    # Quick flow parameter
    a::P = Param(0.9, bounds=(0.1, 10.0), desc="Quickflow coefficient (1/tau_q)")

    # Groundwater parameters (NEW)
    tau_s::P = Param(10.0, bounds=(10.0, 900.0), desc="Slow flow time constant (days)")
    vs::P = Param(0.5, bounds=(0.1, 0.9), desc="Recharge fraction (to GW)")
    L::P = Param(0.0, bounds=(-10.0, 10.0), desc="Constant loss (ML/day)")
    specific_yield::P = Param(0.15, bounds=(0.05, 0.25), desc="Aquifer specific yield")

    # Elevation parameters (fixed, not calibrated)
    gauge_elevation::Float64 = 100.0  # Stream bed elevation (mAHD)
    bore_ground_elevation::Float64 = 110.0  # Bore ground surface (mAHD)

    # State variables
    storage::Array{A} = [100.0]  # CMD (mm)
    quick_store::Array{A} = [0.0]  # Quickflow storage (ML)
    gw_storage::Array{A} = [0.0]  # Groundwater storage (ML)
    effective_rainfall::Array{A} = []  # mm/day
    et::Array{A} = []  # mm/day
    quickflow::Array{A} = []  # ML/day
    baseflow::Array{A} = []  # ML/day
    outflow::Array{A} = []  # ML/day
    inflow::Array{A} = []  # ML/day

    obj_func::Function = obj_func
end


"""
    IHACRESBilinearNodeGW(name::String, spec::AbstractDict)

Construct IHACRESBilinearNodeGW from specification dictionary (YAML/Dict).

# Required spec fields
- `area` : Catchment area (km²)
- `aquifer_area` : Aquifer area (km²) - for depth conversions
- `initial_storage` : Initial CMD (mm)
- `parameters` : Dict with parameter values
- `gauge_elevation` : Stream bed elevation (mAHD)
- `bore_ground_elevation` : Bore ground elevation (mAHD)

# Optional spec fields
- `initial_gw_storage` : Initial GW storage (ML), defaults to 0.0
"""
function IHACRESBilinearNodeGW(name::String, spec::AbstractDict)
    n = create_node(IHACRESBilinearNodeGW, name, spec["area"])
    node_params = copy(spec["parameters"])

    # Set initial CMD storage
    n.storage = [spec["initial_storage"]]

    # Set initial GW storage (if provided)
    if haskey(spec, "initial_gw_storage")
        n.gw_storage = [spec["initial_gw_storage"]]
    end

    # Set aquifer area (required)
    n.aquifer_area = spec["aquifer_area"]

    # Set elevation parameters
    n.gauge_elevation = spec["gauge_elevation"]
    n.bore_ground_elevation = spec["bore_ground_elevation"]

    # Set parameters from spec
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
    IHACRESBilinearNodeGW(name, area, d, d2, e, f, a, tau_s, vs, L, specific_yield,
                          gauge_elev, bore_elev, cmd_store, quick_store, gw_store)

Construct IHACRESBilinearNodeGW with direct parameter specification.
"""
function IHACRESBilinearNodeGW(
    name::String, area::Float64,
    d::Float64, d2::Float64, e::Float64, f::Float64,
    a::Float64, alpha::Float64,
    tau_s::Float64, vs::Float64, L::Float64, specific_yield::Float64,
    gauge_elev::Float64, bore_elev::Float64,
    cmd_store::Float64, quick_store::Float64, gw_store::Float64
)
    n = create_node(IHACRESBilinearNodeGW, name, area)

    # Set parameters
    n.d = Param(d, bounds=n.d.bounds)
    n.d2 = Param(d2, bounds=n.d2.bounds)
    n.e = Param(e, bounds=n.e.bounds)
    n.f = Param(f, bounds=n.f.bounds)
    n.a = Param(a, bounds=n.a.bounds)
    n.alpha = Param(alpha, bounds=n.alpha.bounds)
    n.tau_s = Param(tau_s, bounds=n.tau_s.bounds)
    n.vs = Param(vs, bounds=n.vs.bounds)
    n.L = Param(L, bounds=n.L.bounds)
    n.specific_yield = Param(specific_yield, bounds=n.specific_yield.bounds)

    # Set elevations
    n.gauge_elevation = gauge_elev
    n.bore_ground_elevation = bore_elev

    # Set initial states
    n.storage = [cmd_store]
    n.quick_store = [quick_store]
    n.gw_storage = [gw_store]

    return n
end


"""
    prep_state!(node::IHACRESBilinearNodeGW, timesteps::Int64)::Nothing

Prepare node for simulation by pre-allocating state arrays.
"""
function prep_state!(node::IHACRESBilinearNodeGW, timesteps::Int64)::Nothing
    resize!(node.storage, timesteps + 1)
    node.storage[2:end] .= 0.0

    resize!(node.quick_store, timesteps + 1)
    node.quick_store[2:end] .= 0.0

    resize!(node.gw_storage, timesteps + 1)
    node.gw_storage[2:end] .= 0.0

    node.effective_rainfall = fill(0.0, timesteps)
    node.et = fill(0.0, timesteps)
    node.quickflow = fill(0.0, timesteps)
    node.baseflow = fill(0.0, timesteps)
    node.outflow = fill(0.0, timesteps)
    node.inflow = fill(0.0, timesteps)

    return nothing
end


"""
    reset!(node::IHACRESBilinearNodeGW)::Nothing

Reset node state to initial conditions.
"""
function reset!(node::IHACRESBilinearNodeGW)::Nothing
    node.storage = [node.storage[1]]
    node.quick_store = [node.quick_store[1]]
    node.gw_storage = [node.gw_storage[1]]
    node.effective_rainfall = []
    node.et = []
    node.quickflow = []
    node.baseflow = []
    node.outflow = []
    node.inflow = []

    return nothing
end


"""
    run_timestep!(node::IHACRESBilinearNodeGW, rain::F, evap::F, ts::Int64;
                  inflow=nothing, extraction=nothing, exchange=nothing)::F where {F<:Float64}

Run IHACRESBilinearNodeGW for one timestep.

# Arguments
- `node` : IHACRESBilinearNodeGW instance
- `rain` : Rainfall (mm/day)
- `evap` : Evaporation or temperature
- `ts` : Current timestep index
- `inflow` : Inflow from upstream nodes (ML/day), default 0.0
- `extraction` : Groundwater extraction (ML/day), default 0.0
- `exchange` : External GW flux (ML/day), default 0.0

# Returns
Total outflow (ML/day) = quickflow + baseflow + inflow - extraction
"""
function run_timestep!(
    node::IHACRESBilinearNodeGW,
    rain::F,
    evap::F,
    ts::Int64;
    inflow=nothing,
    extraction=nothing,
    exchange=nothing
)::F where {F<:Float64}
    # Get current states
    current_cmd::F = node.storage[ts]
    quick_store::F = node.quick_store[ts]
    gw_store::F = node.gw_storage[ts]

    # Parse optional inputs
    in_flow::F = timestep_value(ts, node.name, "inflow", inflow)
    gw_ext::F = timestep_value(ts, node.name, "extraction", extraction)
    gw_exchange::F = timestep_value(ts, node.name, "exchange", exchange)

    # 1. Calculate interim CMD and effective rainfall
    (mf, e_rainfall, recharge) = calc_ft_interim_cmd(
        current_cmd,
        rain,
        node.d.val::F,
        node.d2.val::F,
        node.alpha.val::F
    )

    # 2. Calculate ET
    et::F = calc_ET_from_E(
        node.e.val::F,
        evap,
        mf,
        node.f.val::F,
        node.d.val::F
    )

    # 3. Update CMD
    cmd::F = calc_cmd(
        current_cmd,
        rain,
        et,
        e_rainfall,
        recharge
    )

    # 4. Calculate quickflow (uses 1-vs fraction of effective rainfall)
    (nq_store, quickflow) = calc_ft_quick_flow(
        quick_store,
        e_rainfall,
        node.area,
        node.a.val::F,
        node.vs.val::F
    )

    # 5. Calculate GW storage and baseflow (uses vs fraction of effective rainfall)
    (ngw_store, baseflow) = routing(
        gw_store,
        e_rainfall,
        node.area,
        node.tau_s.val::F,
        node.vs.val::F,
        gw_ext,
        node.L.val::F,
        gw_exchange
    )

    # 6. Calculate total outflow
    # Total = local quickflow + local baseflow + upstream inflow - local extraction
    outflow = quickflow + baseflow + in_flow

    # 7. Update state
    node.storage[ts+1] = cmd
    node.quick_store[ts+1] = nq_store
    node.gw_storage[ts+1] = ngw_store

    node.inflow[ts] = in_flow
    node.quickflow[ts] = quickflow
    node.baseflow[ts] = baseflow
    node.outflow[ts] = outflow
    node.effective_rainfall[ts] = e_rainfall
    node.et[ts] = et

    return outflow
end


"""
    param_info(node::IHACRESBilinearNodeGW)::Tuple

Extract parameter names, values, and bounds for calibration.

Returns tuple of (param_names, values, bounds).
"""
function param_info(node::IHACRESBilinearNodeGW)::Tuple
    tmp = Model(node)
    values = collect(tmp[:val])
    bounds = collect(tmp[:bounds])
    param_names = collect(tmp[:fieldname])

    return param_names, values, bounds
end


"""
    update_params!(node::IHACRESBilinearNodeGW, d, d2, e, f, a, alpha, tau_s, vs, L, specific_yield)

Update node parameters (typically after calibration).
"""
function update_params!(
    node::IHACRESBilinearNodeGW,
    d::Float64, d2::Float64, e::Float64, f::Float64,
    a::Float64, alpha::Float64,
    tau_s::Float64, vs::Float64, L::Float64, specific_yield::Float64
)::Nothing
    node.d = Param(d, bounds=node.d.bounds)
    node.d2 = Param(d2, bounds=node.d2.bounds)
    node.e = Param(e, bounds=node.e.bounds)
    node.f = Param(f, bounds=node.f.bounds)
    node.a = Param(a, bounds=node.a.bounds)
    node.alpha = Param(alpha, bounds=node.alpha.bounds)
    node.tau_s = Param(tau_s, bounds=node.tau_s.bounds)
    node.vs = Param(vs, bounds=node.vs.bounds)
    node.L = Param(L, bounds=node.L.bounds)
    node.specific_yield = Param(specific_yield, bounds=node.specific_yield.bounds)

    return nothing
end


"""
    extract_spec!(node::IHACRESBilinearNodeGW, spec::AbstractDict)::Nothing

Extract node-specific details to specification dictionary.
"""
function extract_spec!(node::IHACRESBilinearNodeGW, spec::AbstractDict)::Nothing
    spec["initial_storage"] = node.storage[1]
    spec["initial_gw_storage"] = node.gw_storage[1]
    spec["aquifer_area"] = node.aquifer_area
    spec["gauge_elevation"] = node.gauge_elevation
    spec["bore_ground_elevation"] = node.bore_ground_elevation

    return nothing
end


"""
    get_simulated_depth(node::IHACRESBilinearNodeGW, ts::Int64)::Float64

Get simulated depth below ground surface for comparison with bore observations.

Returns depth in meters below ground surface at bore location.
Positive values indicate water table is below ground surface (normal condition).
"""
function get_simulated_depth(node::IHACRESBilinearNodeGW, ts::Int64)::Float64
    return convert_storage_to_depth(
        node.gw_storage[ts],
        node.aquifer_area,  # Use aquifer_area, not catchment area
        node.specific_yield.val,
        node.gauge_elevation,
        node.bore_ground_elevation
    )
end


"""
    get_water_table_elevation(node::IHACRESBilinearNodeGW, ts::Int64)::Float64

Get water table elevation in mAHD at timestep ts.

Returns absolute elevation of water table in Australian Height Datum (mAHD).
This is the elevation of the water table surface, useful for regional mapping
and comparing water levels between different locations.

# Example
```julia
wt_elevation = get_water_table_elevation(node, 100)  # Get elevation at day 100
println("Water table at ", wt_elevation, " mAHD")
```
"""
function get_water_table_elevation(node::IHACRESBilinearNodeGW, ts::Int64)::Float64
    area_m2 = node.aquifer_area * 1e6  # Use aquifer_area
    water_table_height = (node.gw_storage[ts] * 1000.0) / (area_m2 * node.specific_yield.val)
    return node.gauge_elevation + water_table_height
end


"""
    get_water_table_height_above_streambed(node::IHACRESBilinearNodeGW, ts::Int64)::Float64

Get height of water table above stream bed in meters at timestep ts.

Returns height above gauge elevation (the reference where storage=0).
Positive values indicate water table is above stream bed (normal for gaining stream).
Negative values indicate water table is below stream bed (losing stream or disconnected).

# Example
```julia
height = get_water_table_height_above_streambed(node, 100)
if height > 0
    println("Water table is ", height, " m above stream bed - baseflow likely")
else
    println("Water table is below stream bed - no baseflow")
end
```
"""
function get_water_table_height_above_streambed(node::IHACRESBilinearNodeGW, ts::Int64)::Float64
    area_m2 = node.aquifer_area * 1e6  # Use aquifer_area
    return (node.gw_storage[ts] * 1000.0) / (area_m2 * node.specific_yield.val)
end
