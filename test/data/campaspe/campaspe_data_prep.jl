"""
Script to extract and format data into Streamfall expected format/structure.
"""

using OrderedCollections
using Glob

using Statistics
using CSV, DataFrames, YAML
using Streamfall

# Load climate data - in this case from a CSV file with data for all nodes.
climate_data = CSV.read(
    "climate/climate_historic.csv",
    DataFrame;
    comment="#"
)

# Load in observation data for each gauge in network
outflow_files = readdir(glob"gauges/*_outflow*.csv")
all_flow_data = CSV.read.(outflow_files, DataFrame; comment="#")

# Process data. Later on, flow data will be re-combined into a single dataframe
ordered_flow_data = OrderedDict()
for name in node_names(sn)
    matching_file_pos = findall(occursin.(name, outflow_files))[1]
    suffix = "_outflow_[ML]"

    if name != "406000"
        ordered_flow_data[name] = extract_flow(all_flow_data[matching_file_pos], name, suffix)
        continue
    end

    out = extract_flow(all_flow_data[matching_file_pos], name, suffix)
    ext = all_flow_data[matching_file_pos][:, ["Date", "406000_extractions_[ML]"]]
    df = DataFrame(Date=out.Date)

    # Note that dam releases must have "releases" in their column name.
    df[!, "406000_releases_[ML]"] = out[:, "406000"] .+ ext[:, "406000_extractions_[ML]"]

    ordered_flow_data[name*"_releases_[ML]"] = df
end

# Streamfall provides some convenience functions to align datasets.
# In the next few lines
aligned_climate, aligned_flow... = align_time_frame(climate_data, values(ordered_flow_data)...)
_, extraction_data = align_time_frame(aligned_climate, ordered_flow_data["406000_releases_[ML]"])

# Read in dam levels for calibration
hist_dam_levels = CSV.read(
    "gauges/406000_historic_levels_for_fit.csv",
    DataFrame;
    comment="#"
)
_, aligned_dam_levels = align_time_frame(aligned_climate, hist_dam_levels)

# Combine all data into a single DataFrame
calib_data = hcat(
    aligned_flow[1][:, [:Date]],
    [aligned_flow[i][:, 2:end] for i in 1:length(keys(aligned_flow))]...
)

# For dams, level data should be used for calibration, so we replace that here.
rename!(calib_data, "406000_releases_[ML]"=>"406000")
calib_data[!, "406000"] = aligned_dam_levels[!, "Dam Level [mAHD]"]

# We now have a dataset for calibration (`calib_data`) and a dataset indicating the
# historic dam extractions (`extraction_data`) which is added onto releases to obtain total
# dam outflow.
# `extraction_data` may also hold water extractions at each "reach".

CSV.write("gauges/outflow_and_level.csv", calib_data)
CSV.write("climate/climate.csv", aligned_climate)
CSV.write("gauges/dam_extraction.csv", extraction_data)
