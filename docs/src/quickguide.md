# Quick Start Guide

AC OPF using ADMM

```julia
using PMADA
using Ipopt

model_type = PMADA.ACPPowerModel
run_dopf_admm("matpower/case3.m", model_type, Ipopt.Optimizer; tol=1e-4, max_iteration=1000, verbose = true, alpha=1000)
```

## Getting Results

Results are stored in the data dictionary contains subsystems information
```julia
result = run_dopf_admm("matpower/case3.m", model_type, Ipopt.Optimizer; tol=1e-4, max_iteration=1000, verbose = true, alpha=1000)
```
