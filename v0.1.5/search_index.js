var documenterSearchIndex = {"docs":
[{"location":"metrics/#Available-metrics","page":"Available metrics","title":"Available metrics","text":"","category":"section"},{"location":"metrics/","page":"Available metrics","title":"Available metrics","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"metrics.jl\"]","category":"page"},{"location":"metrics/#Streamfall.ADJ_R2-Tuple{Any, Any, Int64}","page":"Available metrics","title":"Streamfall.ADJ_R2","text":"Determine adjusted R^2\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\np::Int : number of explanatory variables\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.BKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.BKGE","text":"Bounded KGE, bounded between -1 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.BmKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.BmKGE","text":"Bounded modified KGE between -1 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.BnpKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.BnpKGE","text":"Bounded non-parametric KGE between -1 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.KGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.KGE","text":"Calculate the 2009 Kling-Gupta Efficiency (KGE) metric.\n\nA KGE score of 1 means perfect fit. A score < -0.41 indicates that the mean of observations provides better estimates (see Knoben et al., 2019).\n\nNote: Although similar, NSE and KGE cannot be directly compared.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\nReferences\n\nGupta, H.V., Kling, H., Yilmaz, K.K., Martinez, G.F., 2009.  Decomposition of the mean squared error and NSE performance criteria:  Implications for improving hydrological modelling.  Journal of Hydrology 377, 80–91.  https://doi.org/10.1016/j.jhydrol.2009.08.003\nKnoben, W.J.M., Freer, J.E., Woods, R.A., 2019.  Technical note: Inherent benchmark or not? Comparing Nash-Sutcliffe and Kling-Gupta efficiency scores (preprint).  Catchment hydrology/Modelling approaches.  https://doi.org/10.5194/hess-2019-327\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.LME-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.LME","text":"Liu Mean Efficiency metric.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\nReferences\n\nLiu, D., 2020.  A rational performance criterion for hydrological model.  Journal of Hydrology 590, 125488.  https://doi.org/10.1016/j.jhydrol.2020.125488\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.NKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.NKGE","text":"Normalized KGE between 0 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.NNSE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.NNSE","text":"Normalized Nash-Sutcliffe Efficiency score (bounded between 0 and 1).\n\nReferences\n\nNossent, J., Bauwens, W., 2012.  Application of a normalized Nash-Sutcliffe efficiency to improve the accuracy of the Sobol’ sensitivity analysis of a hydrological model.  EGU General Assembly Conference Abstracts 237.\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.NSE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.NSE","text":"The Nash-Sutcliffe Efficiency score\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.NmKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.NmKGE","text":"Normalized modified KGE between 0 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.NnpKGE-Tuple{Array, Array}","page":"Available metrics","title":"Streamfall.NnpKGE","text":"Normalized non-parametric KGE between 0 and 1.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.R2-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.R2","text":"Coefficient of determination (R^2)\n\nAliases NSE()\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.RMSE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.RMSE","text":"Root Mean Square Error\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.mKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.mKGE","text":"Calculate the modified KGE metric (2012).\n\nAlso known as KGE prime (KGE').\n\nArguments\n\nobs::Vector: observations\nsim::Vector : modeled results\n\nReferences\n\nKling, H., Fuchs, M., Paulin, M., 2012.  Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios.  Journal of Hydrology 424–425, 264–277.  https://doi.org/10.1016/j.jhydrol.2012.01.011\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.naive_split_metric-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.naive_split_metric","text":"Naive approach to split metrics.\n\nSplit metrics are a meta-objective optimization approach which segments data into subperiods. The objective function is calculated for each subperiod and then recombined. The approach addresses the lack of consideration of dry years  with least-squares.\n\nIn Fowler et al., [1] the subperiod is one year. This method is \"naive\" in  that the time series is partitioned into N chunks of n_members. Therefore, leap years or partial years are not considered.\n\nArguments\n\nobs::Vector : Historic observations to compare against\nsim::Vector : Modeled time series\nn_members::Int : number of members per chunk, defaults to 365\nmetric::Function : Objective function to apply, defaults to NNSE\ncomb_method::Function : Recombination method, defaults to mean\n\nReferences\n\nFowler, K., Peel, M., Western, A., Zhang, L., 2018.  Improved Rainfall-Runoff Calibration for Drying Climate: Choice of Objective Function.  Water Resources Research 54, 3392–3408.  https://doi.org/10.1029/2017WR022466\n\n\n\n\n\n","category":"method"},{"location":"metrics/#Streamfall.npKGE-Tuple{Any, Any}","page":"Available metrics","title":"Streamfall.npKGE","text":"Calculate the non-parametric Kling-Gupta Efficiency (KGE) metric.\n\nArguments\n\nobs::Vector : observations\nsim::Vector : modeled results\n\nReferences\n\nPool, S., Vis, M., Seibert, J., 2018.  Evaluating model performance: towards a non-parametric variant of the Kling-Gupta efficiency.  Hydrological Sciences Journal 63, 1941–1953.  https://doi.org/10.1080/02626667.2018.1552002\n\n\n\n\n\n","category":"method"},{"location":"simple_showcase/#A-simple-example","page":"A simple example","title":"A simple example","text":"","category":"section"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"In this example we showcase a two-node network to represent water levels at Lake Eppalock,  a dam in the Lower Campaspe catchment located in North-Central Victoria, Australia.","category":"page"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"The map below shows Lake Eppalock, along with relevant gauge locations/data  (click the brown dots on the map to see further gauge details).","category":"page"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"This example uses the setup as detailed in Calibration setup.","category":"page"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"<iframe style=\"width: 720px; height: 600px; border: none;\" src=\"https://nationalmap.gov.au/#share=s-dIbct7mdo25m7ZK2EVr7Koi4cMp\" allowFullScreen mozAllowFullScreen webkitAllowFullScreen></iframe>","category":"page"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"@info \"Running example stream...\"\n\nreset!(sn) # clear any previous runs\n\n# Run the dam node and above\ndam_id, dam_node = get_gauge(sn, \"406000\")\nrun_node!(sn, dam_id, climate; water_order=hist_dam_releases)\n\n# Get performance metrics\nh_data = hist_dam_levels[:, \"Dam Level [mAHD]\"]\nn_data = dam_node.level\n\nrmse_score = Streamfall.RMSE(h_data, n_data)\nnnse_score = Streamfall.NNSE(h_data, n_data)\nnse_score = Streamfall.NSE(h_data, n_data)\n\nrmse = round(rmse_score, digits=4)\nnnse = round(nnse_score, digits=4)\nnse = round(nse_score, digits=4)\n\n@info \"Scores:\" rmse_score nnse_score nse_score\n\n\n# Results of model run\nplot(h_data,\n     legend=:bottomleft,\n     title=\"Calibrated IHACRES\\n(RMSE: $(rmse); NSE: $(nse))\",\n     label=\"Historic\", xlabel=\"Day\", ylabel=\"Dam Level [mAHD]\")\n\nplot!(n_data, label=\"IHACRES\")","category":"page"},{"location":"simple_showcase/","page":"A simple example","title":"A simple example","text":"(Image: )","category":"page"},{"location":"calibration_setup/#Calibration-setup","page":"Calibration setup","title":"Calibration setup","text":"","category":"section"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"The calibration examples all rely on the functions shown here.","category":"page"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"List of metrics provided by Streamfall can be found in Available metrics","category":"page"},{"location":"calibration_setup/#Importing-shared/common-packages","page":"Calibration setup","title":"Importing shared/common packages","text":"","category":"section"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"# Ensure dependent data and packages are available\nusing Statistics, DataFrames, CSV\nusing Distributed, BlackBoxOptim\n\nusing ModelParameters\nusing LightGraphs, MetaGraphs\nusing YAML, Plots\nusing Streamfall","category":"page"},{"location":"calibration_setup/#Load-network-specification","page":"Calibration setup","title":"Load network specification","text":"","category":"section"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"Note that the DATA_PATH is pointing to the test/data/campaspe/ directory.","category":"page"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"# Load and generate stream network\nnetwork = YAML.load_file(joinpath(DATA_PATH, \"campaspe_network.yml\"))\nsn = create_network(\"Example Network\", network)","category":"page"},{"location":"calibration_setup/#Loading-historic-data","page":"Calibration setup","title":"Loading historic data","text":"","category":"section"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"# Load climate data\ndate_format = \"YYYY-mm-dd\"\nclimate_data = DataFrame!(CSV.File(joinpath(data_path, \"climate/climate_historic.csv\"),\n                          comment=\"#\",\n                          dateformat=date_format))\n\ndam_level_fn = joinpath(data_path, \"dam/historic_levels_for_fit.csv\")\ndam_releases_fn = joinpath(data_path, \"dam/historic_releases.csv\")\nhist_dam_levels = DataFrame!(CSV.File(dam_level_fn, dateformat=date_format))\nhist_dam_releases = DataFrame!(CSV.File(dam_releases_fn, dateformat=date_format))\n\n# Subset to same range\nclimate_data, hist_dam_levels, hist_dam_releases = Streamfall.align_time_frame(climate_data, \n                                                                               hist_dam_levels, \n                                                                               hist_dam_releases)\n\n# Create historic data alias\nhist_data = Dict(\n    \"406000\" => hist_dam_levels[:, \"Dam Level [mAHD]\"]\n)\n\n# Create climate object\nclimate = Climate(climate_data, \"_rain\", \"_evap\")","category":"page"},{"location":"calibration_setup/#Example-objective-functions","page":"Calibration setup","title":"Example objective functions","text":"","category":"section"},{"location":"calibration_setup/","page":"Calibration setup","title":"Calibration setup","text":"\"\"\"Calibrate current node.\"\"\"\nfunction obj_func(params, climate, sn, v_id, calib_data::Dict)\n\n    this_node = get_node(sn, v_id)\n    update_params!(this_node, params...)\n\n    # Running next node will run this node\n    Streamfall.run_node!(sn, v_id, climate; water_order=hist_dam_releases)\n\n    n_data = this_node.outflow\n    h_data = calib_data[this_node.node_id]\n\n    # Calculate score (NNSE; 0 to 1)\n    NNSE = Streamfall.NNSE(h_data, n_data)\n\n    # Switch fitness direction as we want to minimize\n    score = 1.0 - NNSE\n\n    # reset to clear stored values\n    reset!(sn)\n\n    return score\nend\n\n\n\"\"\"Example objective function when performance of current node is dependent \non the next node.\n\"\"\"\nfunction obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)\n\n    this_node = get_node(sn, v_id)\n    update_params!(this_node, params...)\n\n    # Run next node which will run this node\n    next_node = get_node(sn, next_vid)\n    releases = calib_data[\"$(next_node.node_id)_releases\"]\n    Streamfall.run_node!(sn, next_vid, climate; water_order=releases)\n\n    # Alias data as necessary\n    if next_node.node_id == \"406000\"\n        n_data = next_node.level\n        h_data = calib_data[next_node.node_id]\n    elseif this_node.node_id == \"406000\"\n        n_data = this_node.level\n        h_data = calib_data[this_node.node_id]\n    else\n        n_data = this_node.outflow\n        h_data = calib_data[this_node.node_id]\n    end\n\n    NNSE = Streamfall.NNSE(h_data, n_data)\n    score = 1.0 - NNSE\n\n    reset!(sn)\n\n    return score\nend\n\n\n\"\"\"Alternative objective function for example. \n\nThis uses a naive split meta-objective function using the Normalized KGE' method.\n\nSee `metrics` page for details.\n\"\"\"\nfunction alt_obj_func(params, climate, sn, v_id, next_vid, calib_data::Dict)\n    this_node = get_node(sn, v_id)\n    update_params!(this_node, params...)\n\n    # Run next node (which will also run this node)\n    Streamfall.run_node!(sn, next_vid, climate; water_order=hist_dam_releases)\n\n    next_node = get_node(sn, next_vid)\n    # Alias data as necessary\n    if next_node.node_id == \"406000\"\n        n_data = next_node.level\n        h_data = calib_data[next_node.node_id]\n    elseif this_node.node_id == \"406000\"\n        n_data = this_node.level\n        h_data = calib_data[this_node.node_id]\n    else\n        n_data = this_node.outflow\n        h_data = calib_data[this_node.node_id]\n    end\n\n    split_NmKGE = Streamfall.naive_split_metric(h_data, n_data; n_members=365, metric=Streamfall.NmKGE, comb_method=mean)\n    score = 1.0 - split_NmKGE\n\n    reset!(sn)\n\n    return score\nend","category":"page"},{"location":"use_methods/#Methods-to-run-a-network-or-node","page":"Methods to run a network or node","title":"Methods to run a network or node","text":"","category":"section"},{"location":"use_methods/","page":"Methods to run a network or node","title":"Methods to run a network or node","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"Streamfall.jl\"]","category":"page"},{"location":"use_methods/#Streamfall.align_time_frame-Union{Tuple{Vararg{T, N} where N}, Tuple{T}} where T<:DataFrames.DataFrame","page":"Methods to run a network or node","title":"Streamfall.align_time_frame","text":"align_time_frame(timeseries::T...)\n\nSubset an arbitrary number of DataFrames to their shared period of time.\n\nReturns subsetted copy of data in same order as input.\n\nExample\n\njulia> climate, streamflow = align_time_frame(climate, streamflow)\n\n\n\n\n\n","category":"method"},{"location":"use_methods/#Streamfall.find_common_timeframe-Union{Tuple{Vararg{T, N} where N}, Tuple{T}} where T<:DataFrames.DataFrame","page":"Methods to run a network or node","title":"Streamfall.find_common_timeframe","text":"find_common_timeframe(timeseries::T...)\n\nFind common time frame between time series.\n\nRequires that all DataFrames have a \"Date\" column.\n\n\n\n\n\n","category":"method"},{"location":"use_methods/#Streamfall.run_catchment!-Tuple{Streamfall.StreamfallNetwork, Climate}","page":"Methods to run a network or node","title":"Streamfall.run_catchment!","text":"run_catchment!(sn::StreamfallNetwork, climate::Climate; water_order=nothing, exchange=nothing)\n\nRun scenario for an entire catchment.\n\n\n\n\n\n","category":"method"},{"location":"use_methods/#Streamfall.run_node!-Tuple{Streamfall.NetworkNode, Any}","page":"Methods to run a network or node","title":"Streamfall.run_node!","text":"run_node!(node::NetworkNode, climate; \n          inflow=nothing, water_order=nothing, exchange=nothing)\n\nRun a specific node, and only that node, for all time steps.\n\nArguments\n\nnode::NetworkNode :\nclimate::Climate :\ninflow::Vector : Time series of inflows from any upstream node.\nwater_order::Vector : Time series of water extractions from this subcatchment\nexchange::Vector : Time series of groundwater flux\n\n\n\n\n\n","category":"method"},{"location":"use_methods/#Streamfall.run_node!-Tuple{Streamfall.StreamfallNetwork, Int64, Climate, Int64}","page":"Methods to run a network or node","title":"Streamfall.run_node!","text":"run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate, timestep::Int; \n          water_order::Union{DataFrame, Nothing}=nothing,\n          exchange::Union{DataFrame, Nothing}=nothing)\n\nRun a model attached to a node for a given time step. Recurses upstream as needed.\n\nArguments\n\nsn::StreamfallNetwork\nnode_id::Int\nclimate::Climate\ntimestep::Int\nwater_order::DataFrame : Volume of water to be extracted (in ML/timestep)\nexchange::DataFrame : Volume of flux (in ML/timestep), where negative values are losses to the groundwater system\n\n\n\n\n\n","category":"method"},{"location":"use_methods/#Streamfall.run_node!-Tuple{Streamfall.StreamfallNetwork, Int64, Climate}","page":"Methods to run a network or node","title":"Streamfall.run_node!","text":"run_node!(sn::StreamfallNetwork, node_id::Int, climate::Climate; \n          water_order=nothing, exchange=nothing)::Nothing\n\nRun model for all time steps, recursing upstream as needed.\n\nArguments\n\nsn::StreamfallNetwork\nnode_id::Int : node to run in the network\nclimate::Climate : Climate object holding rainfall and evaporation data (or temperature)\nwater_order::Vector : water orders for each time step (defaults to nothing)\nexchange::Vector : exchange with groundwater system at each time step (defaults to nothing)\n\n\n\n\n\n","category":"method"},{"location":"nodes/#Implemented-node-types","page":"Implemented node types","title":"Implemented node types","text":"","category":"section"},{"location":"nodes/#IHACRES","page":"Implemented node types","title":"IHACRES","text":"","category":"section"},{"location":"nodes/","page":"Implemented node types","title":"Implemented node types","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"IHACRESNode.jl\"]","category":"page"},{"location":"nodes/#Streamfall.param_info-Tuple{IHACRESNode}","page":"Implemented node types","title":"Streamfall.param_info","text":"param_info(node::IHACRESNode; with_level::Bool = true)::Tuple\n\nExtract node parameter values and bounds\n\n\n\n\n\n","category":"method"},{"location":"nodes/#Streamfall.reset!-Tuple{IHACRESNode}","page":"Implemented node types","title":"Streamfall.reset!","text":"reset!(s_node::IHACRESNode)::Nothing\n\nReset node. Clears all states back to their initial values.\n\n\n\n\n\n","category":"method"},{"location":"nodes/#Streamfall.run_node!","page":"Implemented node types","title":"Streamfall.run_node!","text":"run_node!(s_node::BilinearNode,\n          rain::Float64,\n          evap::Float64,\n          inflow::Float64,\n          ext::Float64,\n          gw_exchange::Float64=0.0;\n          current_store::Union{Nothing, Float64}=nothing,\n          quick_store::Union{Nothing, Float64}=nothing,\n          slow_store::Union{Nothing, Float64}=nothing)::Tuple{Float64, Float64}\n\nRun node with ET data to calculate outflow and update state.\n\nParameters\n\ns_node::BilinearNode : IHACRESNode\nrain : rainfall for time step\nevap : evapotranspiration for time step\ninflow : inflow from previous node\next : irrigation and other water extractions\ngw_exchange : flux in ML where positive is contribution to stream, negative is infiltration\ncurrent_store : replacement cmd value\nquickstore : replacement quickstore value\nslowstore : replacement slowstore value\n\nReturns\n\nfloat, outflow from node, stream level\n\n\n\n\n\n","category":"function"},{"location":"nodes/#Streamfall.run_node_with_temp!","page":"Implemented node types","title":"Streamfall.run_node_with_temp!","text":"run_node_with_temp!(s_node::BilinearNode,\n                    rain::Float64,\n                    temp::Float64,\n                    inflow::Float64,\n                    ext::Float64,\n                    gw_exchange::Float64=0.0;\n                    current_store=nothing,\n                    quick_store=nothing,\n                    slow_store=nothing)::Tuple{Float64, Float64}\n\nRun node with temperature data to calculate outflow and update state.\n\n\n\n\n\n","category":"function"},{"location":"nodes/#Streamfall.update_params!-Tuple{BilinearNode, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}","page":"Implemented node types","title":"Streamfall.update_params!","text":"update_params!(node::BilinearNode, d::Float64, d2::Float64, e::Float64, f::Float64,\n               a::Float64, b::Float64, s_coef::Float64, alpha::Float64)::Nothing\n\nUpdate model parameters.\n\n\n\n\n\n","category":"method"},{"location":"nodes/#Streamfall.update_params!-Tuple{BilinearNode{ModelParameters.Param}, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}","page":"Implemented node types","title":"Streamfall.update_params!","text":"update_params!(node::BilinearNode{Param}, d::Float64, d2::Float64, e::Float64, f::Float64,\n               a::Float64, b::Float64, s_coef::Float64, alpha::Float64,\n               p1::Float64, p2::Float64, p3::Float64, p4::Float64, p5::Float64, p6::Float64, p7::Float64, p8::Float64, CTF::Float64)::Nothing\n\nUpdate all parameters.\n\n\n\n\n\n","category":"method"},{"location":"nodes/#IHACRES-Expuh","page":"Implemented node types","title":"IHACRES - Expuh","text":"","category":"section"},{"location":"nodes/","page":"Implemented node types","title":"Implemented node types","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"IHACRESExpuhNode.jl\"]","category":"page"},{"location":"nodes/#Dam-Level","page":"Implemented node types","title":"Dam Level","text":"","category":"section"},{"location":"nodes/","page":"Implemented node types","title":"Implemented node types","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"DamNode.jl\"]","category":"page"},{"location":"nodes/#Streamfall.param_info-Tuple{DamNode}","page":"Implemented node types","title":"Streamfall.param_info","text":"Extract node parameter values and bounds\n\n\n\n\n\n","category":"method"},{"location":"nodes/#Streamfall.run_node!-2","page":"Implemented node types","title":"Streamfall.run_node!","text":"Calculate outflow for the dam node for a single time step.\n\nParameters\n\nnode : DamNode rain : Float64, rainfall in mm et : Float64, evapotranspiration data in mm irrigext : Float64, irrigation extractions  extractions : Float64, extraction data in ML gwflux : Float64, groundwater interaction\n\n:returns: numeric, outflow from Dam\n\n\n\n\n\n","category":"function"},{"location":"nodes/#Streamfall.update_params!-Tuple{DamNode, Float64}","page":"Implemented node types","title":"Streamfall.update_params!","text":"\n\n\n\n","category":"method"},{"location":"nodes/#Streamfall.update_volume-NTuple{9, Any}","page":"Implemented node types","title":"Streamfall.update_volume","text":"Update dam volume for timestep\n\nParameters\n\nvolume : float, current water volume in ML nodeinflow : float, inflow from previous node in ML gamma : float, groundwater exchange (positive is gain from gw flow, negative is loss to infiltration) rain : float, rainfall input evap : float, evaporation loss infiltration : float, infiltration loss area : float, dam surface area in square kilometers extractions : float, water extraction from dam in ML discharge : float, discharge from dam in ML maxstore : float, maximum dam storage in ML\n\nReturns\n\nfloat, volume of water stored in dam\n\n\n\n\n\n","category":"method"},{"location":"calibration/#Example-calibration","page":"Example calibration","title":"Example calibration","text":"","category":"section"},{"location":"calibration/","page":"Example calibration","title":"Example calibration","text":"# Import common packages and functions\n# This is shown under Calibration Setup\ninclude(\"_obj_func_definition.jl\")\n\n\n\"\"\"Example calibration function.\n\nIllustrate model calibration using the BlackBoxOptim package.\n\"\"\"\nfunction calibrate(sn, v_id, climate, calib_data)\n\n    ins = inlets(sn, v_id)\n\n    # Recurse through and calibrate all nodes upstream\n    if !isempty(ins)\n        for nid in ins\n            calibrate(sn, nid, climate, calib_data)\n        end\n    end\n\n    this_node = get_node(sn, v_id)\n\n    # Create new optimization function (see definition inside Calibration Setup)\n    opt_func = x -> obj_func(x, climate, sn, v_id, calib_data)\n\n    # Get node parameters (default values and bounds)\n    x0, param_bounds = param_info(this_node; with_level=false)\n    opt = bbsetup(opt_func; SearchRange=param_bounds,\n                  Method=:adaptive_de_rand_1_bin_radiuslimited,\n                  MaxTime=300.0,  # time in seconds to spend\n                  TraceInterval=30.0,\n                  PopulationSize=100,\n                  # Workers=workers()\n                  )\n    \n    res = bboptimize(opt)\n\n    bs = best_candidate(res)\n    @info \"Calibrated $(v_id) ($(this_node.node_id)), with score: $(best_fitness(res))\"\n    @info \"Best Params:\" collect(bs)\n\n    # Update node with calibrated parameters\n    update_params!(this_node, bs...)\n\n    return res, opt\nend\n\n\nv_id, node = get_gauge(sn, \"406219\")\n@info \"Starting calibration...\"\nres, opt = calibrate(sn, v_id, climate, hist_data)\n\nbest_params = best_candidate(res)\n\n@info best_fitness(res)\n@info best_params\n\n\nupdate_params!(node, best_params...)\nStreamfall.run_node!(sn, v_id, climate)\n\nh_data = hist_data[\"406219\"]\nn_data = node.outflow\n\n@info \"Outflow NNSE:\" Streamfall.NNSE(h_data, n_data)\n@info \"Outflow RMSE:\" Streamfall.RMSE(h_data, n_data)\nreset!(node)\n\ndam_id, dam_node = get_gauge(sn, \"406000\")\ntimesteps = sim_length(climate)\nfor ts in (1:timesteps)\n    run_node!(sn, dam_id, climate, ts; water_order=hist_dam_releases)\nend\n\nh_data = hist_dam_levels[:, \"Dam Level [mAHD]\"]\nn_data = dam_node.level\n\nrmse = Streamfall.RMSE(h_data, n_data)\nnse = Streamfall.NSE(h_data, n_data)\n\n# Results of model run\nplot(h_data,\n     legend=:bottomleft,\n     title=\"Calibrated IHACRES\\n(RMSE: $(rmse); NSE: $(nse))\",\n     label=\"Historic\", xlabel=\"Day\", ylabel=\"Dam Level [mAHD]\")\n\nplot!(n_data, label=\"IHACRES\")\n\nsavefig(\"calibrated_example.png\")\n\n# 1:1 Plot\nscatter(h_data, n_data, legend=false, \n        markerstrokewidth=0, markerstrokealpha=0, alpha=0.2)\nplot!(h_data, h_data, color=:red, markersize=.1, markerstrokewidth=0,\n      xlabel=\"Historic [mAHD]\", ylabel=\"IHACRES [mAHD]\", title=\"Historic vs Modelled\")\n\nsavefig(\"calibration_1to1.png\")\n\n\n# NNSE: 0.9643; RMSE: 1.43553\n# d: 84.28015146853407\n# d2: 2.4224106535469145\n# e: 0.8129590022893607\n# f: 2.579276454391652\n# a: 5.923379062122229\n# b: 0.0989925603647026\n# storage_coef: 1.8613364808233752  # gw storage factor\n# alpha: 0.7279050097363565","category":"page"},{"location":"calibration/","page":"Example calibration","title":"Example calibration","text":"(Image: )","category":"page"},{"location":"calibration/","page":"Example calibration","title":"Example calibration","text":"(Image: )","category":"page"},{"location":"#Streamfall.jl-Documentation","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"","category":"section"},{"location":"","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"Streamfall: An experimental graph-based streamflow modelling system written in Julialang.","category":"page"},{"location":"","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"Aims of the project are to leverage the Julia language and ecosystem to allow/enable:","category":"page"},{"location":"","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"Quick application and exploratory analysis\nUse of different rainfall-runoff models in tandem [aspiration]\nModelling and assessment of interacting socio-environmental systems\nParallel scenario runs","category":"page"},{"location":"","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"Note: the only model currently available is the IHACRES rainfall-runoff model, leveraging ihacres_nim.","category":"page"},{"location":"","page":"Streamfall.jl Documentation","title":"Streamfall.jl Documentation","text":"LightGraphs and MetaGraphs are used underneath for network traversal/analysis.","category":"page"},{"location":"multisystem_showcase/#A-simple-showcase-of-multi-system-considerations","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"","category":"section"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"This page is a draft.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Here we showcase a two-node network representing a river and a dam downstream.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"The Lower Campaspe catchment - a small semi-arid basin in North-Central Victoria, Australia - is used for the example here.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Figure of catchment","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"As a graph, the network looks like this:","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Figure of two-node network","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"The dam is the primary water store for farmers in the area but is also used for recreational activities (camping, boating, fishing, etc) by local enthusiasts and vacationers. The Campaspe river is also home to a culturally and ecologically significant population of fish. A certain level of flow must be ensured at key times during the year to support and maintain their population levels.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"In this hypothetical study, local stakeholders would like to have an idea of the range of possible dam levels under a range of environmental watering policies, farmer water use, and how this may impact the level of enjoyment by vacationers.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"The possible environmental watering strategies are defined as:","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Implicit watering:   No purposeful releases for environmental demands. Assume natural inflows and agricultural water orders provide sufficient water flow for ecological purposes.\nExplicit watering:   Assume agricultural water orders partially fulfill environmental needs. Water is released as needed to meet any deficit.\nPrioritized watering:   Water for environmental purposes are prioritised and are managed separately from agricultural demands.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"For the purpose of this example, the farm water requirements are given as a volume of daily water releases throughout a growing season; the period of time over which a crop can grow. This figure may be provided by another model in practice. The growing season is assumed to be between X and Y, with the daily water requirements over that period being between X and Y.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"An index value is used to provide indications of the suitability of dam levels for recreational purposes.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"explain how recreational index works","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Another indicator model is used to show how often environmental needs are met.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"explain how the environmental indicator model works","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"An overview of the system under investigation can then be conceptualized like:","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Conceptual figure of the system","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Where water flows into the dam, and water is released to fulfill water needs of the users downstream. Note that \"water users\" as defined here includes the environment itself.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"First, we define a two-node Streamfall Network which represents the river and dam:","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"We can then generate a number of scenarios representing a mix of the management strategies and water demands, as listed above.","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"# Code to generate scenarios","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"# code showing how to run the model(s)","category":"page"},{"location":"multisystem_showcase/","page":"A simple showcase of multi-system considerations","title":"A simple showcase of multi-system considerations","text":"Analysis and wrap up...","category":"page"},{"location":"network/#Defining-a-network","page":"Defining a network","title":"Defining a network","text":"","category":"section"},{"location":"network/","page":"Defining a network","title":"Defining a network","text":"Modules = [Streamfall]\nOrder   = [:function, :type]\nPages   = [\"Network.jl\"]","category":"page"},{"location":"network/#Streamfall.create_network-Tuple{String, Dict}","page":"Defining a network","title":"Streamfall.create_network","text":"create_network(name::String, network::Dict)::StreamfallNetwork\n\nCreate a StreamNetwork from a YAML-derived specification.\n\nExample\n\njulia> network_spec = YAML.load_file(\"example_network.yml\")\njulia> sn = create_network(\"Example Network\", network_spec)\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.create_node-Tuple{Streamfall.StreamfallNetwork, String, Dict, Int64}","page":"Defining a network","title":"Streamfall.create_node","text":"create_node(sn::StreamfallNetwork, node_name::String, details::Dict, nid::Int)\ncreate_node(mg::MetaGraph, node_name::String, details::Dict, nid::Int)\n\nCreate a node specified by with given name (if it does not exist).\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.find_inlets_and_outlets-Tuple{Streamfall.StreamfallNetwork}","page":"Defining a network","title":"Streamfall.find_inlets_and_outlets","text":"Find all inlets and outlets in a network.\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.in_or_out-Tuple{Any, Any}","page":"Defining a network","title":"Streamfall.in_or_out","text":"Determine a node's connection\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.inlets-Tuple{Streamfall.StreamfallNetwork, String}","page":"Defining a network","title":"Streamfall.inlets","text":"inlets(sn::StreamfallNetwork, node_id::String)\n\nFind nodes which provides inflows for given node.\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.outlets-Tuple{Streamfall.StreamfallNetwork, String}","page":"Defining a network","title":"Streamfall.outlets","text":"outlets(sn::StreamfallNetwork, node_id::String)\n\nFind node immediately downstream from given node.\n\n\n\n\n\n","category":"method"},{"location":"network/#Streamfall.reset!-Tuple{Streamfall.StreamfallNetwork}","page":"Defining a network","title":"Streamfall.reset!","text":"reset!(sn::StreamfallNetwork)::Nothing\n\nReset a network.\n\n\n\n\n\n","category":"method"}]
}