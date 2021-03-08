using YAML
using Streamfall


network = YAML.load_file("../test/data/campaspe/campaspe_network.yml")
g, mg = create_network("Example Network", network)

inlets, outlets = find_inlets_and_outlets(g)

@info "Network has the following inlets and outlets:" inlets outlets
