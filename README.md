# Thermodynamic ranking of pathways in reaction networks

This github repository contains the code used to generate the results presented in the paper: [Thermodynamic ranking of pathways in reaction networks](https://doi.org/10.48550/arXiv.2506.23496)

## Dependencies

The code is written in the programming language *[julia](https://julialang.org/)* (version 1.10.4).
The following packages are needed to generate and output the data:

- [JuMP](https://jump.dev/JuMP.jl/stable/)
- [Ipopt](https://github.com/jump-dev/Ipopt.jl)
- DelimitedFiles (included in default *julia* installation)

Additionally the following packages are needed for the plots:

- [Plots](https://docs.juliaplots.org/stable/)
- StatsPlots (included in the *Plots* package)
- Colors (included in the *Plots* package)
- [LaTeXStrings](https://github.com/JuliaStrings/LaTeXStrings.jl)

## Overview

The subfolders of this repository contain the following scripts:

- *4species*: contains the scripts related to the 4-species system
- *5species*: contains the scripts related to the 4-species system
- *ac_compete*: contains the scripts related to the competing autocatalytic-cycle system

