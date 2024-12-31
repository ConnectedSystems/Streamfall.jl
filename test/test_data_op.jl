using Statistics
using Test

using Dates, DataFrames, CSV
using Streamfall

@testset "Data alignment" begin

    here = @__DIR__
    climate_data = joinpath(here, "data/campaspe/climate/climate_historic.csv")
    dam_data_loc = joinpath(here, "data/campaspe/gauges")

    # Read in test data
    climate_data = CSV.File(climate_data, comment="#", dateformat="YYYY-mm-dd") |> DataFrame
    hist_dam_levels = CSV.File(joinpath(dam_data_loc, "406000_historic_levels_for_fit.csv"), dateformat="YYYY-mm-dd") |> DataFrame
    hist_dam_releases = CSV.File(joinpath(dam_data_loc, "406000_historic_outflow.csv"), dateformat="YYYY-mm-dd") |> DataFrame

    @testset "Min/max date finding" begin
        min_date, max_date = Streamfall.find_common_timeframe(climate_data, hist_dam_levels, hist_dam_releases)

        @test min_date == Date("1981-01-01")
        @test max_date == Date("2015-07-30")
    end

    @testset "Aligning mismatched datasets" begin
        # Mismatched start/end dates
        dam_levels = hist_dam_levels[10:end, :]
        dam_releases = hist_dam_releases[20:10000, :]

        # Subset to common time period
        c, dl, dr = Streamfall.align_time_frame(climate_data, dam_levels, dam_releases)

        # Ensure time periods are the same
        @test c.Date[1] == dl.Date[1] == dr.Date[1]
        @test c.Date[end] == dl.Date[end] == dr.Date[end]

        # Ensure originals have not been modified
        @test climate_data.Date[1] != c.Date[1]
        @test hist_dam_levels.Date[1] != dl.Date[1]
        @test hist_dam_releases.Date[1] != dr.Date[1]

        @test hist_dam_levels.Date[end] != dl.Date[end]
        @test hist_dam_releases.Date[end] != dr.Date[end]
    end
end
