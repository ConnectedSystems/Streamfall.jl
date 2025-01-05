module GraphMakieExt

using Streamfall
import Streamfall: StreamfallNetwork

using Makie

using Graphs
using GraphMakie


"""
    plot_network(sn::StreamfallNetwork)

Simple plot of stream network.
"""
function Streamfall.NetworkViz.plot_network(sn::StreamfallNetwork)
    node_labels = ["$(sn[i].name)\n"*string(nameof(typeof(sn[i]))) for i in vertices(sn.mg)]
    f, ax, sp = graphplot(
        sn.mg;
        nlabels=node_labels,
        nlabels_align=(:center, :bottom),
        node_size=16,
        node_color=:blue
    )
    # f.scene.padding = (50, 50, 50, 50)  # (left, right, bottom, top)

    # Turn off grid lines
    ax.xgridvisible = false
    ax.ygridvisible = false

    # Hide ticks and labels
    hidedecorations!(ax)

    reset_limits!(ax)

    return f
end

end