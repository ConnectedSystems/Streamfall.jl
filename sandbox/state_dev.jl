using Streamfall
using Statistics, CSV, DataFrames
using Distributed, BlackBoxOptim
using Base.Iterators
using Plots
import Streamfall.Analysis: TemporalCrossSection


include("IHACRESCMDStateNode.jl")


function v_metric(obs, sim; sp=365, metric=Streamfall.NmKGE)
    res = Streamfall.naive_split_metric(obs, sim, sp, metric)
    m = mean(res)
    cv = std(res) / m

    V = (m*0.9) + ((ℯ^-cv)*0.1)  #  (ℯ .- 1.5)

    return V
end


function v_mod_metric(obs, sim; sp=365, metric=Streamfall.NmKGE)
    res = Streamfall.naive_split_metric(obs, sim, sp, metric)
    m = mean(res)
    sd = std(res)

    V = (m*0.9) + ((ℯ^-sd)*0.1)

    return V
end


# load observations
HERE = @__DIR__
DATA_PATH = joinpath(HERE, "../test/data/cotter/")

date_format = "YYYY-mm-dd"
obs_data = CSV.File(joinpath(DATA_PATH, "climate/CAMELS-AUS_410730.csv"),
                    comment="#",
                    dateformat=date_format) |> DataFrame

hist_streamflow = obs_data[:, ["Date", "410730_Q"]]
hist_streamflow = rename(hist_streamflow, "410730_Q"=>"410730")

climate_data = obs_data[:, ["Date", "410730_P", "410730_PET"]]
climate = Climate(climate_data, "_P", "_PET")

state_node = create_node(CMDStateNode, "410730", 129.2)


c_date = climate_data.Date
obs = hist_streamflow[:, "410730"]
figs = []
outflows = []
display_timeframe = 13330:(13330+(365*1))


burn_in = 1826

split_mod_v = (obs, sim) -> 1.0 - v_mod_metric(obs[burn_in:end], sim[burn_in:end]; sp=14, metric=Streamfall.NnpKGE)


# calibrate!(state_node, climate, hist_streamflow; metric=split_mod_v, MaxTime=10.0, TraceMode=:silent)

res, opt = calibrate!(state_node, climate, hist_streamflow; metric=split_mod_v, MaxTime=10.0, TraceMode=:silent);

