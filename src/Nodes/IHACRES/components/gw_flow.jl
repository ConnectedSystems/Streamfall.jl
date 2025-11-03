"""
    routing(
        prev_gw_storage::Float64,
        recharge::Float64,
        area::Float64,
        tau_s::Float64,
        vs::Float64,
        extraction::Float64,
        loss::Float64,
        gw_exchange::Float64
    )::Tuple{Float64,Float64}

Calculate groundwater storage and baseflow using IHACRES_GW formulation.

Based on Ivkovic et al. (2005) groundwater module extension to IHACRES.

# Arguments
- `prev_gw_storage` : Groundwater storage at t-1 (ML)
- `recharge` : Recharge from effective rainfall (mm/day)
- `area` : Catchment area (km²)
- `tau_s` : Slow flow time constant (days), controls recession rate
- `vs` : Recharge fraction (0-1), proportion of effective rainfall to GW
- `extraction` : Groundwater extraction/pumping (ML/day)
- `loss` : Constant loss term L (ML/day), positive = loss, negative = gain
- `gw_exchange` : External GW flux (ML/day), positive = inflow, negative = outflow

# Returns
Tuple of (gw_storage, baseflow) both in ML/day

# References
- Ivkovic et al. (2005a) - IHACRES_GW model description
- Croke & Jakeman (2004) - IHACRES-CMD module

# Notes
- When gw_storage ≤ 0, baseflow = 0 (water table below stream bed)
- The coefficient 'a' is derived from tau_s: a = 1/τ_s (exponential recession coefficient)
- Recharge is converted from mm/day to ML/day by: recharge_ML = vs * recharge * area
- Storage and baseflow are partitioned as: G[t] = G_tmp/(1+a), Q[t] = a*G_tmp/(1+a)
"""
function routing(
    prev_gw_storage::Float64,
    recharge::Float64,
    area::Float64,
    tau_s::Float64,
    vs::Float64,
    extraction::Float64,
    loss::Float64,
    gw_exchange::Float64
)::Tuple{Float64,Float64}
    # Calculate baseflow coefficient 'a' from time constant tau_s
    # For exponential recession: a = 1/τ_s
    # This gives Q = a*G which produces exponential decay
    a = 1.0 / tau_s

    # Convert recharge from mm/day to ML/day
    # vs is the fraction that goes to groundwater (remainder goes to quick flow)
    recharge_ML = vs * recharge * area

    # Update storage (before baseflow calculation)
    # Storage increases with recharge and gw_exchange (inflow)
    # Storage decreases with extraction and loss
    gw_tmp = prev_gw_storage + recharge_ML - extraction - loss + gw_exchange

    # Calculate baseflow and final storage
    # Baseflow only occurs if storage is positive (water table above stream bed)
    if gw_tmp > 0.0
        # Split storage into baseflow (discharge) and remaining storage
        # From IHACRES_GW: Q_s = a * G_t and G_t = G_tmp / (1 + a)
        baseflow = a * gw_tmp / (1.0 + a)
        gw_storage = gw_tmp / (1.0 + a)
    else
        # Water table below stream bed - no baseflow
        baseflow = 0.0
        gw_storage = max(0.0, gw_tmp)  # Prevent negative storage
    end

    # Ensure baseflow is non-negative
    baseflow = max(0.0, baseflow)

    return (gw_storage, baseflow)
end


"""
    calc_ft_quick_flow(
        prev_quick::Float64,
        e_rain::Float64,
        area::Float64,
        a::Float64,
        vs::Float64
    )::Tuple{Float64,Float64}

Calculate quickflow store and discharge.

Extracted from calc_ft_flows to allow separate quick and slow flow calculation
in IHACRES_GW where slow flow is replaced by groundwater module.

# Arguments
- `prev_quick` : Previous quickflow storage (ML)
- `e_rain` : Effective rainfall (mm/day)
- `area` : Catchment area (km²)
- `a` : Quickflow storage coefficient (1/τ_q)
- `vs` : Recharge fraction (proportion going to GW, remainder goes to quick)

# Returns
Tuple of (quick_store, quickflow) in ML/day

# Notes
- Only (1 - vs) fraction of effective rainfall goes to quickflow
- vs fraction goes to groundwater recharge instead
"""
function calc_ft_quick_flow(
    prev_quick::Float64,
    e_rain::Float64,
    area::Float64,
    a::Float64,
    vs::Float64
)::Tuple{Float64,Float64}
    # Partition effective rainfall between quickflow and GW recharge
    # (1 - vs) goes to quick flow, vs goes to GW recharge
    quick_rain = (1.0 - vs) * e_rain

    # Calculate quickflow store and discharge
    tmp_calc = max(0.0, prev_quick + (quick_rain * area))

    if tmp_calc > 0.0
        alpha = exp(-a)
        beta = (1.0 - alpha) * tmp_calc
        quick_store = alpha * tmp_calc
        quickflow = beta
    else
        quick_store = tmp_calc
        quickflow = 0.0
    end

    @assert quickflow >= 0.0 "Quickflow cannot be negative: $quickflow"

    return (quick_store, quickflow)
end


"""
    convert_storage_to_depth(
        gw_storage_ML::Float64,
        area_km2::Float64,
        specific_yield::Float64,
        stream_bed_elevation_mAHD::Float64,
        bore_ground_elevation_mAHD::Float64
    )::Float64

Convert groundwater storage (ML) to depth below ground surface (m).

Used for comparing simulated water table with bore observations of "Depth to water (m)".

# Arguments
- `gw_storage_ML` : Groundwater storage volume (ML)
- `area_km2` : Catchment area (km²)
- `specific_yield` : Aquifer specific yield (dimensionless, typically 0.1-0.3)
- `stream_bed_elevation_mAHD` : Stream bed elevation (mAHD), reference where storage=0
- `bore_ground_elevation_mAHD` : Ground surface elevation at bore (mAHD)

# Returns
Depth below ground surface (m), positive values indicate water table below ground

# Calculation Steps
1. Convert storage to water table height above stream bed:
   height = storage (ML) * 1000 / (area_m² * specific_yield)
2. Calculate water table elevation:
   elevation = stream_bed + height
3. Calculate depth below surface:
   depth = ground_surface - water_table_elevation

# Notes
- When storage = 0: water table is at stream bed level
- Positive depth = water table below ground (normal condition)
- Negative depth = water table above ground (flooding/artesian)
"""
function convert_storage_to_depth(
    gw_storage_ML::Float64,
    area_km2::Float64,
    specific_yield::Float64,
    stream_bed_elevation_mAHD::Float64,
    bore_ground_elevation_mAHD::Float64
)::Float64
    # Convert area from km² to m²
    area_m2 = area_km2 * 1e6

    # Convert storage (ML) to water table height (m) above stream bed
    # Storage in ML → m³ (×1000) → divide by (area × specific_yield)
    water_table_height_m = (gw_storage_ML * 1000.0) / (area_m2 * specific_yield)

    # Calculate water table elevation in mAHD
    # When storage = 0, water table is at stream bed elevation
    water_table_elevation_mAHD = stream_bed_elevation_mAHD + water_table_height_m

    # Calculate depth below ground surface (positive = below ground)
    depth_below_surface_m = bore_ground_elevation_mAHD - water_table_elevation_mAHD

    return depth_below_surface_m
end


"""
    convert_depth_to_storage(
        depth_below_surface_m::Float64,
        area_km2::Float64,
        specific_yield::Float64,
        stream_bed_elevation_mAHD::Float64,
        bore_ground_elevation_mAHD::Float64
    )::Float64

Convert depth below ground surface (m) to groundwater storage (ML).

Inverse of convert_storage_to_depth. Useful for initializing storage from observed bore depths.

# Arguments
- `depth_below_surface_m` : Observed depth to water below ground (m)
- `area_km2` : Catchment area (km²)
- `specific_yield` : Aquifer specific yield (dimensionless)
- `stream_bed_elevation_mAHD` : Stream bed elevation (mAHD)
- `bore_ground_elevation_mAHD` : Ground surface elevation at bore (mAHD)

# Returns
Groundwater storage (ML)

# Notes
- Used to back-calculate initial storage from first bore observation
- Can return negative storage if water table is below stream bed
"""
function convert_depth_to_storage(
    depth_below_surface_m::Float64,
    area_km2::Float64,
    specific_yield::Float64,
    stream_bed_elevation_mAHD::Float64,
    bore_ground_elevation_mAHD::Float64
)::Float64
    # Calculate water table elevation from depth
    water_table_elevation_mAHD = bore_ground_elevation_mAHD - depth_below_surface_m

    # Calculate height above stream bed
    water_table_height_m = water_table_elevation_mAHD - stream_bed_elevation_mAHD

    # Convert to storage
    area_m2 = area_km2 * 1e6
    gw_storage_ML = (water_table_height_m * area_m2 * specific_yield) / 1000.0

    return gw_storage_ML
end
