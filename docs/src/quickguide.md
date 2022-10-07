# Quick Start Guide

AC OPF using ADMM

```julia
using PowerModelsADA
using Ipopt

model_type = ACPPowerModel
run_dopf_admm("test/data/case_RTS.m", model_type, Ipopt.Optimizer; verbose=1)
```

## Getting Results

Results are stored in the data dictionary contains subsystems information
```julia
result = run_dopf_admm("test/data/case_RTS.m", model_type, Ipopt.Optimizer; verbose=1)
```
