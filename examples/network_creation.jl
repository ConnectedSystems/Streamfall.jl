using YAML
using Streamfall

data_dir = joinpath(dirname(dirname(pathof(Streamfall))), "test/data")
network = YAML.load_file(joinpath(data_dir, "campaspe/campaspe_network.yml"))
sn = create_network("Example Network", network)

inlets, outlets = find_inlets_and_outlets(sn)

@info "Network has the following inlets and outlets:" inlets outlets
