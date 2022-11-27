#  PowerModelsADA.jl 

Status:
[![CI](https://github.com/mkhraijah/PowerModelsADA.jl/workflows/CI/badge.svg)](https://github.com/mkhraijah/PowerModelsADA.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/mkhraijah/PowerModelsADA.jl/branch/main/graph/badge.svg?token=371LK4OBZG)](https://codecov.io/gh/mkhraijah/PowerModelsADA.jl)
[![Documentation](https://github.com/mkhraijah/PowerModelsADA.jl/workflows/Documentation/badge.svg)](https://mkhraijah.github.io/PowerModelsADA.jl/)
</p>

## Overview

[`PowerModelsADA.jl`](https://github.com/mkhraijah/PowerModelsADA.jl) (Power Models Alternating Distributed Algorithms) provides a framework to solve Optimal Power Flow (OPF) problems using alternating distributed algorithms. The package allows to use different distributed algorithms. `PowerModelsADA` is built on top of [`PowerModels.jl`](https://github.com/lanl-ansi/PowerModels.jl) and [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) to model and solve the subproblems.


## Distributed Algorithms 
The `PowerModelsADA` framework is designed to easily incorporate new alternating distributed algorithms. The framework provides means to decompose a test case into multiple areas, model the subproblems associated with each area using `PowerModels`, solve the supropblems in parallel using multi-threading, communicate the shared data between the areas, and calculate the mismatches to decide if the termination criteria are satisfied.


The current version of `PowerModelsADA` implements four distributed algorithms: 

- Alternating Direction Method of Multipliers (ADMM)
- Analytical Target Cascading (ATC)
- Auxiliary Problem Principle (APP)
- Augmented Lagrangian Alternating Direction Inexact Newton (ALADIN)

`PowerModelsADA` can be extended to include variations of the existing algorithms or new user-defined algorithms. 

<!--
 More details about the formulations and algorithm implementations are shown in [Technical Specifications](https://mkhraijah.github.io/PowerModelsADA.jl/dev/specification/)

## Installation

`PowerModelsADA` can be installed using the Julia package manager with

```julia
using Pkg
Pkg.add("PowerModelsADA")
```  
-->



## Examples

An example demonstrating how to code up and solve the OPF problem with distributed algorithms is found in [Quick Start Guide](https://mkhraijah.github.io/PowerModelsADA.jl/dev/quickguide/) section of the documentation.
