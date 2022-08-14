# PMADA.jl

```@meta
CurrentModule = PMADA
```
## Overview

[PMADA.jl](https://github.com/mkhraijah/PMADA.jl) (Power Models Alternating Distributed Algorithms) provides a framework to solve the Optimal Power Flow (OPF) problem using alternating distributed algorithms. The package allows to use different distributed algorithms such as Alternating Direction Method of Multipliers (ADMM) or user-defined algorithms. PMADA is built on top of `PowerModels.jl` to define and solve the subproblems.

## Installation

PMADA can be installed using the Julia package manager with

```julia
] add PMADA
```

## Examples

An example demonstrating how to code up and solve the OPF problem with distributed algorithms is found in [Quick Start Guide](@ref quickguide.md) section of the documentation.
