using Streamfall


function test_run(s::Float64, e::Float64)::Array{Float64, 1}

    a = Float64[0.0, 0.0]
    @ccall ihacres.test_ptr(a::Ptr{Cdouble}, s::Cdouble, e::Cdouble)::Cvoid

    return a
end

@info test_run(100.0, 200.0)

test_node = IHACRESNode{Float64}(
    100.0,  # area
    100.0,  # d 
    0.1,  # d2
    0.1,  # e 
    0.1,  # f
    54.35,  # a
    0.012,  # b
    2.5,  # storage_coef
    0.240195,  # alpha
    [100.0],  # storage
    [100.0],  # quickflow
    [100.0],  # slowflow
    [],  # outflow
    [],  # effective_rainfall
    [],  # et
    []   # inflow
)


@info "Outflow:" run_node(test_node, 6.0, 3.0, 50.0, 10.0)

@info test_node.slowflow