# PowerModelsADA.jl

```@meta
CurrentModule = PowerModelsADA
```

## Overview

[`PowerModelsADA.jl`](https://github.com/mkhraijah/PowerModelsADA.jl) (Power Models Alternating Distributed Algorithms) provides a framework to solve Optimal Power Flow (OPF) problems using alternating distributed algorithms. The package allows to use different distributed algorithms. `PowerModelsADA` is built on top of [`PowerModels.jl`](https://github.com/lanl-ansi/PowerModels.jl) and [`JuMP.jl`](https://github.com/jump-dev/JuMP.jl) to model and solve the subproblems.

## Distributed Algorithms

The `PowerModelsADA` framework is designed to easily incorporate new alternating distributed algorithms. The framework provides means to decompose a test case into multiple areas, model the subproblems associated with each area using `PowerModels`, solve the supropblems in parallel using multi-threading or multi-processing via [`Distributed Computing`](https://docs.julialang.org/en/v1/manual/distributed-computing/), communicate the shared data between the areas, and calculate the mismatches to decide if the termination criteria are satisfied.

The current version of `PowerModelsADA` implements four distributed algorithms:

- Alternating Direction Method of Multipliers (ADMM)
- Analytical Target Cascading (ATC)
- Auxiliary Problem Principle (APP)
- Augmented Lagrangian Alternating Direction Inexact Newton (ALADIN)

`PowerModelsADA` can be extended to include variations of the existing algorithms or new user-defined algorithms. More details about the formulations and algorithm implementations are shown in [Technical Specifications](https://mkhraijah.github.io/PowerModelsADA.jl/dev/specification/)

## Installation

`PowerModelsADA` can be installed using the Julia package manager with

```julia
using Pkg
Pkg.add("PowerModelsADA")
```

## Examples

An example demonstrating how to code up and solve the OPF problem with distributed algorithms is found in [Quick Start Guide](https://mkhraijah.github.io/PowerModelsADA.jl/dev/quickguide/) section of the documentation.

## Contributions

Contributions and enhancements of `PowerModelADA` are welcomed and encouraged. Please feel free to fork this repository and share your contributions to the main branch with a pull request.

## Citation

If you find `PowerModelsADA` useful for your work, please cite our [paper](https://ieeexplore.ieee.org/document/10262198):

```bibtex
@ARTICLE{alkhraijah2023powermodelsada,
  author={Alkhraijah, Mohannad and Harris, Rachel and Coffrin, Carleton and Molzahn, Daniel K.},
  journal={IEEE Transactions on Power Systems}, 
  title={PowerModelsADA: A Framework for Solving Optimal Power Flow using Distributed Algorithms}, 
  year={2023},
  volume={},
  number={},
  pages={1-4},
  doi={10.1109/TPWRS.2023.3318858}
}
```

## Acknowledgments

This work is partially supported by the NSF AI Institute for Advances in Optimization (Award #2112533).
