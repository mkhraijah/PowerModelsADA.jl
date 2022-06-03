
path = joinpath(@__DIR__,"..")
cd(path)
using Pkg
Pkg.activate(".")
using Ipopt

##

include("../src/DistributedPowerModels.jl")
DPM = DistributedPowerModels

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
data = DPM.parse_file(case_path)


## Settings
alpha = 1000
max_iteration = 1000
tol = 1e-4
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)


##  Distributed algorithm settings
model_type = DPM._PM.DCPPowerModel

data_area = DPM.run_dopf_admm(data, model_type, optimizer, alpha=1000, tol=1e-4, max_iteration=1000, verbose = true)

data_area = DPM.run_dopf_app(data, model_type, optimizer, alpha=1000, tol=1e-4, max_iteration=1000, verbose = true)

data_area = DPM.run_dopf_atc(data, model_type, optimizer, alpha=1.05, tol=1e-4, max_iteration=1000, verbose = true)
