using DataFrames

abstract type NetworkNode{A} end

@def network_node begin
    node_id::String
    area::Float64  # area in km^2
end


# function get_climate_data(node::NetworkNode, climate_data::DataFrame)
#     tgt::String = node.node_id
#     rain_prefix = "pr_"
#     et_prefix = "wvap_"

#     rain = climate_data[Symbol(rain_prefix*tgt)]
#     et = climate_data[Symbol(et_prefix*tgt)]

#     return rain, et
# end
