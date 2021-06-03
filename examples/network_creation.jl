using YAML
using Streamfall


network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
sn = create_network("Example Network", network)

inlets, outlets = find_inlets_and_outlets(sn)

@info "Network has the following inlets and outlets:" inlets outlets
