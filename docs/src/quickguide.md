# Quick Start Guide

To solve the OPF problem using the ADMM

```julia
using PowerModelsADA
using Ipopt

model_type = ACPPowerModel
run_dopf_admm("test/data/case_RTS.m", model_type, Ipopt.Optimizer; verbose=1)
```

## Getting Results

The solve function stores the result in a data dictionary contains subsystems information
```julia
result = run_dopf_admm("test/data/case_RTS.m", model_type, Ipopt.Optimizer; print_level=1)
```
