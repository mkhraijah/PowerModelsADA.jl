
path = joinpath(@__DIR__,"..")
using Pkg
Pkg.activate(".")
using Ipopt
using SeDuMi
## Call Packages
include("DPM.jl")

## Read case function with partition file return dictionary of the case
case_path = "data/case300.m"
data = DPM.parse_file(case_path)

partition_path= "data/case300_3areas.csv"
DPM.assign_area!(data, partition_path)

## Settings
alpha = 100;
max_iteration = 10000;
tol = 1e-2;
distributed_algorithm = "APP" # Select algorithm among APP ADMM ATC
optimizer = Gurobi.Optimizer
pf_model = "DC" # Select model among DC AC SOC SDP QC

settings = DPM.set_setting(distributed_algorithm, tol, max_iteration)

## Run distirbuted opf
DPM.run_dopf(data, pf_model, optimizer, settings, alpha, verbose = true)
