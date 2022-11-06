using Plots, StatsPlots


using DataFrames, CSV, Serialization
using Streamfall, BlackBoxOptim, Bootstrap, Statistics
import Dates: Date, month, monthday, week, dayofyear


HERE = @__DIR__
DATA_PATH = joinpath(HERE, "../test/data/cotter/")

# # Load observations
# date_format = "YYYY-mm-dd"
# obs_data = CSV.File(joinpath(DATA_PATH, "climate/CAMELS-AUS_410730.csv"),
#                         comment="#",
#                         dateformat=date_format) |> DataFrame

# hist_streamflow = obs_data[:, ["Date", "410730_Q"]]
# climate_data = obs_data[:, ["Date", "410730_P", "410730_PET"]]
# climate = Climate(climate_data, "_P", "_PET")

# burn_in = 366  # 1 year burn-in period + 1 day
# obs = hist_streamflow[:, "410730_Q"]
# log_obs = log.(obs .+ 0.02)

valid_range = 14354:22704


prb_obs = open("problem_obs.txt", "r") do f
    x = readlines(f)[1]
    x = split.(x, ",")
    parse.(Float64, string.(x))
end

prb_sim = open("problem_sim.txt", "r") do f
    x = readlines(f)[1]
    x = split.(x, ",")
    parse.(Float64, string.(x))
end

prb_dt = open("problem_dates.txt", "r") do f
    x = readlines(f)[1]
    x = split.(x, ",")
    Date.(x)
end

replace!(prb_obs, 0=>1e-6)
replace!(prb_sim, 0=>1e-6)

uc_xsect = Streamfall.temporal_cross_section(prb_dt, prb_obs, prb_sim; show_extremes=false, yscale=:log10)  # 

display(uc_xsect)
