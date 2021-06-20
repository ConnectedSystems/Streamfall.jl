# Streamfall.jl Documentation

Streamfall: An experimental graph-based streamflow modelling system written in [Julialang](http://julialang.org/).

Aims of the project are to leverage the Julia language and ecosystem to allow/enable:
- Quick application and exploratory analysis
- Use of different rainfall-runoff models in tandem [**aspiration**]
- Modelling and assessment of interacting socio-environmental systems
- Parallel scenario runs

**Note:** the only model currently available is the IHACRES rainfall-runoff model, leveraging [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim).

[LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used underneath for network traversal/analysis.
