using YAML
using Test
using Streamfall, MetaGraphs


@testset "Bare node creation" begin
    test_node = BilinearNode{Float64}(;
        node_id="Test",
        area=100.0,
        route=false
    )

    @info "Outflow:" run_node!(test_node, 6.0, 3.0, 50.0, 10.0)
end


@testset "NaN outputs" begin
    test_node = BilinearNode(
        "Test",  # name/id
        1985.73,  # area
        false,    # route
        200.0,  # d
        2.0,  # d2
        1.0,  # e
        1.675,  # f
        54.35254,  # a
        0.187,  # b
        2.9,  # storage_coef
        0.727,  # alpha
        100.0,  # storage
        0.0,  # quickflow
        0.0  # slowflow
    )

    rain = 7.96848605e+01
    evap = 3.32467909e+00
    inflow = 9.19583373e-11
    extraction = 1.59313987e-12
    gw_exchange = 2.47076926e-11
    loss = 5.51086184e-11
    current_store = 1.03467364e+02
    quick_store = 8.46269687e+02
    slow_store = 3.67133471e+02

    res = run_node!(test_node, rain, evap, inflow, extraction, gw_exchange, loss; 
                    current_store=current_store, quick_store=quick_store, slow_store=slow_store)

    @test !any(isnan, res)
    @info res
end



@testset "Network creation" begin
    # Ensure specified parameter values are being assigned on node creation
    # Load and generate stream network
    network = YAML.load_file("data/campaspe/campaspe_network.yml")
    mg, g = create_network("Example Network", network)

    target_node = get_prop(mg, 1, :node)

    @test target_node.area == 1985.73
    @test target_node.level_params[1] == -3.3502
end


@testset "Interim CMD" begin
    params = (214.6561105573191, 76.6251447, 200.0, 2.0, 0.727)
    current_store, rain, d, d2, alpha = params

    interim_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_interim_cmd(interim_results::Ptr{Cdouble},
                                       current_store::Cdouble,
                                       rain::Cdouble,
                                       d::Cdouble,
                                       d2::Cdouble,
                                       alpha::Cdouble)::Cvoid
    
    (mf, e_rainfall, recharge) = interim_results

    @test !isnan(mf)
    @test !isnan(e_rainfall)
    @test !isnan(recharge)
end

@testset "Catchment Moisture Deficit" begin
    cmd = 100.0
    et = 6.22
    e_rain = 6.83380027058404E-06
    recharge = 3.84930005080411E-06
    rain = 0.0000188

    n_cmd = @ccall IHACRES.calc_cmd(cmd::Cdouble, rain::Cdouble, et::Cdouble, e_rain::Cdouble, recharge::Cdouble)::Float64

    @test isapprox(n_cmd, 106.22, atol=0.001)
end


@testset "IHACRES calculations" begin
    area = 1985.73
    a = 54.352
    b = 0.187
    e_rain = 3.421537294474909e-6
    recharge = 3.2121031313153022e-6
    loss = 0.0

    prev_quick = 100.0
    prev_slow = 100.0

    flow_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_flows(
        flow_results::Ptr{Cdouble},
        prev_quick::Cdouble,
        prev_slow::Cdouble,
        e_rain::Cdouble,
        recharge::Cdouble,
        area::Cdouble,
        a::Cdouble,
        b::Cdouble,
        loss::Cdouble
    )::Cvoid

    @info flow_results

    @test flow_results[1] == (1.0 / (1.0 + a) * (prev_quick + (e_rain * area)))

    b2 = 1.0
    slow_store = prev_slow + (recharge * area) - (loss * b2)
    slow_store = 1.0 / (1.0 + b) * slow_store
    @test flow_results[2] == slow_store

    e_rain = 0.0
    recharge = 0.0

    prev_quick = 3.3317177943791187
    prev_slow = 144.32012122323678

    flow_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_flows(
        flow_results::Ptr{Cdouble},
        prev_quick::Cdouble,
        prev_slow::Cdouble,
        e_rain::Cdouble,
        recharge::Cdouble,
        area::Cdouble,
        a::Cdouble,
        b::Cdouble,
        loss::Cdouble
    )::Cvoid

    @info flow_results

    @test flow_results[1] == (1.0 / (1.0 + a) * (prev_quick + (e_rain * area)))

    b2 = 1.0
    slow_store = prev_slow + (recharge * area) - (loss * b2)
    slow_store = 1.0 / (1.0 + b) * slow_store
    @test flow_results[2] == slow_store
end