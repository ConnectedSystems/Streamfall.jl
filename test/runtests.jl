using Test
using Streamfall

@testset "Bare node creation" begin
    test_node = IHACRESNode{Float64}(;
        node_id="Test",
        area=100.0
    )

    @info "Outflow:" run_node!(test_node, 6.0, 3.0, 50.0, 10.0)
end
