using YAML
using Test
using Streamfall, MetaGraphs

@testset "Bare node creation" begin
    test_node = IHACRESNode{Float64}(;
        node_id="Test",
        area=100.0
    )

    @info "Outflow:" run_node!(test_node, 6.0, 3.0, 50.0, 10.0)
end


@testset "Network creation" begin
    # Ensure specified parameter values are being assigned on node creation
    # Load and generate stream network
    network = YAML.load_file("data/campaspe/campaspe_network.yml")
    g, mg = create_network("Example Network", network)

    target_node = get_prop(mg, 1, :node)

    @test target_node.a == 54.352

    @test target_node.level_params[1] == -3.3502
end


@testset "Interim CMD" begin
    params = (214.6561105573191, 76.6251447, 200.0, 2.0, 0.727)
    current_store, rain, d, d2, alpha = params

    interim_results = [0.0, 0.0, 0.0]
    @ccall IHACRES.calc_ft_interim(interim_results::Ptr{Cdouble},
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
        flow_results::Ref{Cdouble},
        prev_quick::Cdouble,
        prev_slow::Cdouble,
        e_rain::Cdouble,
        recharge::Cdouble,
        area::Cdouble,
        a::Cdouble,
        b::Cdouble,
        loss::Cdouble
    )::Cvoid

    @test flow_results[1] == (1.0 / (1.0 + a) * (100.0 + (e_rain * area)))
end