# Streamfall.jl Documentation

Streamfall: A graph-based streamflow modelling system written in [Julialang](http://julialang.org/).

Streamfall leverages the Julia language and ecosystem to provide:
- Quick hetrogenous modelling of a stream network
- Use of different rainfall-runoff models and their ensembles in tandem
- Modelling and assessment of interacting systems
- A wide range of performance metrics

This package includes implementations of the following:
- GR4J
- HyMod
- IHACRES
- SYMHYD

Performance is expected to be similar to implementations in C and Fortran.

 <table align="center">
  <tr>
    <th>Model</th>
    <th>Full name</th>
    <th>Reference</th>
  </tr>
  <tr>
    <td>GR4J</td>
    <td>modèle du Génie Rural à 4 paramètres Journalier</td>
    <td>
    Perrin, C., Michel, C., Andréassian, V., 2003.
    Improvement of a parsimonious model for streamflow simulation.
    Journal of Hydrology 279, 275-289.
    https://doi.org/10.1016/S0022-1694(03)00225-7
    </td>
  </tr>
  <tr>
    <td>HyMod</td>
    <td>HYdrological MODel</td>
    <td>
    Wagener, T., Boyle, D. P., Lees, M. J., Wheater, H. S., Gupta, H. V.,
    and Sorooshian, S.: A framework for development and applica-
    tion of hydrological models, Hydrol. Earth Syst. Sci., 5, 13–26,
    doi:10.5194/hess-5-13-2001, 2001.
    </td>
  </tr>
  <tr>
    <td>IHACRES</td>
    <td>Identification of unit Hydrographs And Component flows from Rainfall, Evaporation and Streamflow</td>
    <td>
    Croke, B.F.W., Jakeman, A.J. 2004
    A catchment moisture deficit module for the IHACRES rainfall-runoff model,
    Environmental Modelling & Software, 19(1), pp. 1–5.
    doi: 10.1016/j.envsoft.2003.09.001
    </td>
  </tr>
  <tr>
    <td>SYMHYD</td>
    <td>-</td>
    <td>
    Chiew, F. H. S., Peel, M. C., Western, A. W., Singh, V. P., & Frevert, D. (2002). Application and testing of the simple rainfall-runoff model SIMHYD. Mathematical models of small watershed hydrology and applications, 335-367.
    </td>
  </tr>
</table>