module MakieExt

using Statistics
using DataFrames
using Dates

using Streamfall
import Streamfall: TemporalCrossSection

using Makie

function Streamfall.NetworkViz.plot(sn::StreamfallNetwork)
    return graphplot(mg)
end

end