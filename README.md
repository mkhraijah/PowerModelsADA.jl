#  PMADA.jl 

Status:
[![CI](https://github.com/mkhraijah/PMADA.jl/workflows/CI/badge.svg)](https://github.com/mkhraijah/PMADA.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/mkhraijah/PMADA.jl/branch/main/graph/badge.svg?token=371LK4OBZG)](https://codecov.io/gh/mkhraijah/PMADA.jl)
[![Documentation](https://github.com/mkhraijah/PMADA.jl/workflows/Documentation/badge.svg)](https://mkhraijah.github.io/PMADA.jl/)
</p>


## Overview

[`PMADA.jl`](https://github.com/mkhraijah/PMADA.jl) (Power Models Alternating Distributed Algorithms) provides a framework to solve the Optimal Power Flow (OPF) problem using alternating distributed algorithms. The package allows to use different distributed algorithms such as Alternating Direction Method of Multipliers (ADMM) or user-defined algorithms. `PMADA.jl` is built on top of [`PowerModels.jl`](https://github.com/lanl-ansi/PowerModels.jl) to define and solve the subproblems.

## Installation

`PMADA.jl` can be installed using the Julia package manager as follow:

```julia
using Pkg
Pkg.add("PMADA")
```

## Examples

An example demonstrating how to code up and solve the OPF problem with distributed algorithms is found in [Quick Start Guide](https://mkhraijah.github.io/PMADA.jl/dev/quickguide/) section of the documentation.
