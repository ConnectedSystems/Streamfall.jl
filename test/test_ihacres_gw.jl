using Test
using Streamfall


@testset "IHACRES_GW Core Functions" begin

    @testset "routing" begin
        # Test basic GW storage and baseflow calculation

        # Parameters
        prev_storage = 1000.0  # ML
        recharge = 5.0  # mm/day
        area = 254.0  # km²
        tau_s = 20.0  # days
        vs = 0.5  # fraction to GW
        extraction = 10.0  # ML/day
        loss = 0.0  # ML/day
        gw_exchange = 0.0  # ML/day

        (new_storage, baseflow) = routing(
            prev_storage, recharge, area, tau_s, vs, extraction, loss, gw_exchange
        )

        # Check that storage and baseflow are non-negative
        @test new_storage >= 0.0
        @test baseflow >= 0.0

        # Check water balance (approximately)
        recharge_ML = vs * recharge * area
        # Storage change should approximately equal: inflows - outflows
        # new_storage ≈ prev_storage + recharge_ML - baseflow - extraction
        # (approximately, because baseflow depends on storage)
        storage_change = new_storage - prev_storage
        net_flux = recharge_ML - baseflow - extraction
        @test abs(storage_change - net_flux) < 1.0  # Within 1 ML tolerance

        # Test zero storage case (no baseflow)
        (zero_storage, zero_baseflow) = routing(
            0.0, 0.0, area, tau_s, vs, 0.0, 0.0, 0.0
        )
        @test zero_storage == 0.0
        @test zero_baseflow == 0.0

        # Test negative storage input (should be clamped to zero)
        (neg_storage, neg_baseflow) = routing(
            -100.0, 0.0, area, tau_s, vs, 0.0, 0.0, 0.0
        )
        @test neg_storage == 0.0  # Storage clamped to zero (cannot be negative)
        @test neg_baseflow == 0.0  # No baseflow when storage is zero
    end


    @testset "calc_ft_quick_flow" begin
        # Test quickflow calculation

        prev_quick = 100.0  # ML
        e_rain = 10.0  # mm/day
        area = 254.0  # km²
        a = 0.9  # quickflow coefficient
        vs = 0.5  # fraction to GW (so 0.5 to quick)

        (new_quick, quickflow) = calc_ft_quick_flow(prev_quick, e_rain, area, a, vs)

        # Check non-negative
        @test new_quick >= 0.0
        @test quickflow >= 0.0

        # Check that quickflow uses (1-vs) fraction
        expected_quick_rain = (1.0 - vs) * e_rain
        # new_quick should be based on expected_quick_rain
        @test new_quick > 0.0  # Should have some storage

        # Test zero case
        (zero_quick, zero_qflow) = calc_ft_quick_flow(0.0, 0.0, area, a, vs)
        @test zero_quick == 0.0
        @test zero_qflow == 0.0
    end


    @testset "convert_storage_to_depth" begin
        # Test storage to depth conversion

        storage = 5000.0  # ML
        area = 254.0  # km²
        sy = 0.15  # specific yield
        stream_bed = 120.0  # mAHD
        bore_ground = 123.0  # mAHD

        depth = convert_storage_to_depth(storage, area, sy, stream_bed, bore_ground)

        # Depth should be positive (water table below ground)
        @test depth > 0.0

        # Manual calculation check
        area_m2 = area * 1e6
        wt_height = (storage * 1000.0) / (area_m2 * sy)
        wt_elev = stream_bed + wt_height
        expected_depth = bore_ground - wt_elev
        @test abs(depth - expected_depth) < 1e-6

        # Test zero storage (water table at stream bed)
        depth_zero = convert_storage_to_depth(0.0, area, sy, stream_bed, bore_ground)
        @test abs(depth_zero - (bore_ground - stream_bed)) < 1e-6

        # Test negative depth (water table above ground - unlikely but possible)
        large_storage = 50000.0  # ML
        depth_neg = convert_storage_to_depth(large_storage, area, sy, stream_bed, bore_ground)
        # Could be negative if storage is very large
    end


    @testset "convert_depth_to_storage" begin
        # Test depth to storage conversion (inverse operation)

        depth = 10.0  # m below ground
        area = 254.0  # km²
        sy = 0.15
        stream_bed = 120.0  # mAHD
        bore_ground = 123.0  # mAHD

        storage = convert_depth_to_storage(depth, area, sy, stream_bed, bore_ground)

        # Storage could be negative if depth is large (water table below stream)
        # But should be calculable
        @test !isnan(storage)

        # Test round-trip conversion
        storage_test = 5000.0  # ML
        depth_from_storage = convert_storage_to_depth(
            storage_test, area, sy, stream_bed, bore_ground
        )
        storage_back = convert_depth_to_storage(
            depth_from_storage, area, sy, stream_bed, bore_ground
        )
        @test abs(storage_back - storage_test) < 1e-3  # Within 1e-3 ML
    end

end


@testset "IHACRESBilinearNodeGW Node" begin

    @testset "Node construction from dict" begin
        # Test creating node from specification dictionary

        spec = Dict(
            "area" => 254.0,
            "aquifer_area" => 254.0,
            "initial_storage" => 100.0,
            "initial_gw_storage" => 1000.0,
            "gauge_elevation" => 120.0,
            "bore_ground_elevation" => 123.0,
            "parameters" => Dict(
                "d" => 200.0,
                "d2" => 2.0,
                "e" => 1.0,
                "f" => 0.8,
                "a" => 0.9,
                "alpha" => 0.1,
                "tau_s" => 20.0,
                "vs" => 0.5,
                "L" => 0.0,
                "specific_yield" => 0.15
            )
        )

        node = IHACRESBilinearNodeGW("test_node", spec)

        @test node.name == "test_node"
        @test node.area == 254.0
        @test node.storage[1] == 100.0
        @test node.gw_storage[1] == 1000.0
        @test node.gauge_elevation == 120.0
        @test node.bore_ground_elevation == 123.0
        @test node.tau_s.val == 20.0
        @test node.vs.val == 0.5
    end


    @testset "Node construction direct" begin
        # Test creating node with direct parameters

        node = IHACRESBilinearNodeGW(
            "test_direct",
            254.0,  # area
            200.0, 2.0, 1.0, 0.8,  # d, d2, e, f
            0.9, 0.1,  # a, alpha
            20.0, 0.5, 0.0, 0.15,  # tau_s, vs, L, sy
            120.0, 123.0,  # gauge_elev, bore_elev
            100.0, 0.0, 1000.0  # cmd, quick, gw stores
        )

        @test node.name == "test_direct"
        @test node.area == 254.0
        @test node.tau_s.val == 20.0
    end


    @testset "prep_state!" begin
        # Test state preparation

        node = IHACRESBilinearNodeGW(
            "test_prep",
            254.0,
            200.0, 2.0, 1.0, 0.8,
            0.9, 0.1,
            20.0, 0.5, 0.0, 0.15,
            120.0, 123.0,
            100.0, 0.0, 1000.0
        )

        timesteps = 365
        prep_state!(node, timesteps)

        @test length(node.storage) == timesteps + 1
        @test length(node.gw_storage) == timesteps + 1
        @test length(node.quick_store) == timesteps + 1
        @test length(node.outflow) == timesteps
        @test length(node.baseflow) == timesteps
        @test length(node.quickflow) == timesteps
    end


    @testset "reset!" begin
        # Test reset functionality

        node = IHACRESBilinearNodeGW(
            "test_reset",
            254.0,
            200.0, 2.0, 1.0, 0.8,
            0.9, 0.1,
            20.0, 0.5, 0.0, 0.15,
            120.0, 123.0,
            100.0, 0.0, 1000.0
        )

        prep_state!(node, 10)

        # Run a few timesteps to modify state
        node.storage[2] = 90.0
        node.gw_storage[2] = 1100.0

        reset!(node)

        @test length(node.storage) == 1
        @test length(node.gw_storage) == 1
        @test node.storage[1] == 100.0  # Back to initial
        @test node.gw_storage[1] == 1000.0  # Back to initial
    end


    @testset "run_timestep!" begin
        # Test running a single timestep

        node = IHACRESBilinearNodeGW(
            "test_run",
            254.0,
            200.0, 2.0, 1.0, 0.8,
            0.9, 0.1,
            20.0, 0.5, 0.0, 0.15,
            120.0, 123.0,
            100.0, 0.0, 1000.0
        )

        prep_state!(node, 1)

        rain = 10.0  # mm
        evap = 5.0  # mm
        ts = 1

        outflow = run_timestep!(node, rain, evap, ts)

        # Check that outflow was calculated
        @test outflow >= 0.0
        @test node.outflow[1] == outflow

        # Check that states were updated
        @test node.storage[2] != node.storage[1]
        @test node.quickflow[1] >= 0.0
        @test node.baseflow[1] >= 0.0

        # Check that quickflow + baseflow contributes to outflow
        @test node.outflow[1] >= node.quickflow[1] + node.baseflow[1]
    end


    @testset "get_simulated_depth" begin
        # Test getting simulated depth for calibration

        node = IHACRESBilinearNodeGW(
            "test_depth",
            254.0,
            200.0, 2.0, 1.0, 0.8,
            0.9, 0.1,
            20.0, 0.5, 0.0, 0.15,
            120.0, 123.0,
            100.0, 0.0, 5000.0  # 5000 ML GW storage
        )

        depth = get_simulated_depth(node, 1)

        @test depth > 0.0  # Should be positive (below ground)
        @test !isnan(depth)

        # Check consistency with convert_storage_to_depth
        expected_depth = convert_storage_to_depth(
            node.gw_storage[1],
            node.area,
            node.specific_yield.val,
            node.gauge_elevation,
            node.bore_ground_elevation
        )
        @test abs(depth - expected_depth) < 1e-6
    end

end


@testset "Parameter extraction" begin

    node = IHACRESBilinearNodeGW(
        "test_params",
        254.0,
        200.0, 2.0, 1.0, 0.8,
        0.9, 0.1,
        20.0, 0.5, 0.0, 0.15,
        120.0, 123.0,
        100.0, 0.0, 1000.0
    )

    param_names, values, bounds = param_info(node)

    @test length(param_names) == length(values)
    @test length(param_names) == length(bounds)
    @test length(param_names) == 10  # 10 calibratable parameters

    # Check that key parameters are present
    @test :tau_s in param_names
    @test :vs in param_names
    @test :L in param_names
    @test :specific_yield in param_names

end


@testset "Water Balance Partitioning" begin
    # Test that effective rainfall is correctly partitioned between quickflow and GW
    # This is critical for ensuring the model conserves water properly

    @testset "Effective rainfall partitioning" begin
        node = IHACRESBilinearNodeGW(
            "test_balance",
            100.0,  # area
            50.0, 2.0, 1.0, 0.8,  # d, d2, e, f (low d to generate effective rainfall)
            0.5, 0.5,  # a, alpha
            20.0, 0.6, 0.0, 0.15,  # tau_s, vs=0.6, L=0, sy
            120.0, 123.0,
            10.0, 0.0, 1000.0  # low initial CMD to generate effective rainfall
        )

        prep_state!(node, 10)

        # Run a few timesteps with significant rainfall
        rain_events = [50.0, 30.0, 40.0, 20.0, 10.0]
        evap = 5.0

        for (i, rain) in enumerate(rain_events)
            run_timestep!(node, rain, evap, i)
        end

        # Check water balance for days with effective rainfall
        vs = node.vs.val
        area = node.area

        for t in 1:length(rain_events)
            eff_rain = node.effective_rainfall[t]

            if eff_rain > 0.01  # Only check days with significant effective rainfall
                # Convert effective rainfall to ML
                eff_rain_ML = eff_rain * area

                # Check quickflow
                # Note: quickflow includes previous storage, so we need to account for that
                # Expected addition to quick from this timestep: (1-vs) * eff_rain * area
                expected_quick_rain_ML = (1.0 - vs) * eff_rain * area

                # Check GW recharge
                # Expected addition to GW from this timestep: vs * eff_rain * area
                # But the code uses: vs * recharge * area, where recharge != eff_rain
                expected_gw_rain_ML = vs * eff_rain * area

                # Calculate actual GW inflow (storage change + baseflow + L + extraction)
                storage_change = node.gw_storage[t+1] - node.gw_storage[t]
                baseflow = node.baseflow[t]
                L = node.L.val
                actual_gw_inflow = storage_change + baseflow + L  # No extraction in this test

                # Print diagnostic info
                println("\nTimestep $t:")
                println("  Effective rainfall: $(round(eff_rain, digits=2)) mm = $(round(eff_rain_ML, digits=2)) ML")
                println("  Expected to quick ($(round(1-vs, digits=2))*eff_rain): $(round(expected_quick_rain_ML, digits=2)) ML")
                println("  Expected to GW ($(round(vs, digits=2))*eff_rain): $(round(expected_gw_rain_ML, digits=2)) ML")
                println("  Actual to GW (storage change + baseflow + L): $(round(actual_gw_inflow, digits=2)) ML")
                println("  Difference (actual - expected): $(round(actual_gw_inflow - expected_gw_rain_ML, digits=2)) ML")
                println("  Ratio (actual/expected): $(round(actual_gw_inflow / expected_gw_rain_ML, digits=3))")

                # TEST: The partitioning should satisfy vs * eff_rain going to GW
                # Allow 5% tolerance for numerical issues
                @test isapprox(actual_gw_inflow, expected_gw_rain_ML, rtol=0.05)
            end
        end
    end

    @testset "Total water conservation" begin
        # Test that water is conserved over entire simulation
        node = IHACRESBilinearNodeGW(
            "test_conservation",
            100.0,
            50.0, 2.0, 1.0, 0.8,
            0.5, 0.5,
            20.0, 0.5, 0.0, 0.15,  # vs=0.5 for equal split
            120.0, 123.0,
            10.0, 0.0, 1000.0
        )

        prep_state!(node, 20)

        # Run simulation
        for t in 1:20
            rain = 20.0 + 10.0 * sin(t / 3.0)  # Variable rainfall
            evap = 5.0
            run_timestep!(node, rain, evap, t)
        end

        # Sum all water fluxes
        total_eff_rain = sum(node.effective_rainfall) * node.area
        total_quickflow = sum(node.quickflow)
        total_baseflow = sum(node.baseflow)

        # Storage changes
        gw_storage_change = node.gw_storage[end] - node.gw_storage[1]
        quick_storage_change = node.quick_store[end] - node.quick_store[1]

        # Total water out + storage change should equal total effective rainfall
        total_out = total_quickflow + total_baseflow
        total_storage_change = gw_storage_change + quick_storage_change

        println("\nWater Balance Summary:")
        println("  Total effective rainfall: $(round(total_eff_rain, digits=1)) ML")
        println("  Total quickflow: $(round(total_quickflow, digits=1)) ML")
        println("  Total baseflow: $(round(total_baseflow, digits=1)) ML")
        println("  GW storage change: $(round(gw_storage_change, digits=1)) ML")
        println("  Quick storage change: $(round(quick_storage_change, digits=1)) ML")
        println("  Total out + storage change: $(round(total_out + total_storage_change, digits=1)) ML")
        println("  Difference: $(round(total_eff_rain - (total_out + total_storage_change), digits=1)) ML")
        println("  Percentage: $(round(100*(total_out + total_storage_change)/total_eff_rain, digits=1))%")

        # Water balance check (allow 1% error)
        @test isapprox(total_eff_rain, total_out + total_storage_change, rtol=0.01)
    end
end
