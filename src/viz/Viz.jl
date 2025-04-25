module Viz


function plot end
function plot! end
function quickplot end

function temporal_cross_section end

function plot_residuals end

function save_figure end
function save_figure! end


"""Symmetrical log values.

https://kar.kent.ac.uk/32810/2/2012_Bi-symmetric-log-transformation_v5.pdf
https://discourse.julialang.org/t/symmetrical-log-plot/45709/3
"""
function symlog(y)
    return sign.(y) .* log10.(1.0 .+ abs.(y))
end

export plot, plot!, quickplot
export temporal_cross_section, plot_residuals
export save_figure, save_figure!


end