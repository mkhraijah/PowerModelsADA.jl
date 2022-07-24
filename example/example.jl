
using Pkg
Pkg.activate(".")
using Ipopt
using JuMP
## Import PMADA package

using PMADA

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
data = PMADA.parse_file(case_path)


## Settings and
max_iteration = 1000
tol = 1e-4
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

##  Parameter and Power Flow Model selection
alpha = 1000
model_type = PMADA.DCPPowerModel


##  Distributed algorithm settings

data_area = PMADA.run_dopf_admm(data, model_type, optimizer, tol=1e-4, max_iteration=1000, verbose = true, alpha=1000)

data_area = PMADA.run_dopf_app(data, model_type, optimizer, tol=1e-4, max_iteration=1000, verbose = true, alpha=1000)

data_area = PMADA.run_dopf_atc(data, model_type, optimizer, tol=1e-4, max_iteration=1000, verbose = true, alpha=1.05)
