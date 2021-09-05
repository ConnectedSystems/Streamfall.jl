# Streamfall.jl Documentation

Streamfall: An experimental graph-based streamflow modelling system written in [Julialang](http://julialang.org/).

Aims of the project are to leverage the Julia language and ecosystem to allow/enable:
- Quick application and exploratory analysis
- Use of different rainfall-runoff models in tandem
- Modelling and assessment of interacting socio-environmental systems [**aspiration**]
- Parallel scenario runs

Streamfall currently includes interfaces to:
* IHACRES, leveraging [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim)
* HyMOD
* GR4J

[LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl) and [MetaGraphs](https://github.com/JuliaGraphs/MetaGraphs.jl) are used underneath for network traversal/analysis.
