module GraphPlotExt

using Plots, GraphPlot
using Streamfall
import Streamfall: vertices
import Streamfall: StreamfallNetwork

"""
    plot_network(sn::StreamfallNetwork)

Simple plot of stream network.
"""
function Streamfall.NetworkViz.plot_network(sn::StreamfallNetwork; as_html::Bool=false)
    node_labels = ["$(sn[i].name)\n" * string(nameof(typeof(sn[i]))) for i in vertices(sn.mg)]

    if as_html
        plot_func = gplothtml
    else
        plot_func = gplot
    end

    plot_func(sn.mg, nodelabel=node_labels)
end

"""
    save_figure(sn::StreamfallNetwork, fn::String)

Save a figure of the network in SVG format.
"""
function Streamfall.NetworkViz.save_figure(sn::StreamfallNetwork, fn::String)
    draw(SVG(fn, 16cm, 16cm), plot_network(sn))
end

end
