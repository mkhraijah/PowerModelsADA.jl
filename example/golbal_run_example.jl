
path = joinpath(@__DIR__,"..")
using Pkg
Pkg.activate(".")
using Ipopt

## Call Packages
include("../src/DistributedPowerModels.jl")
DPM = DistributedPowerModels

## Read case function with partition file return dictionary of the case
case_path = "test/data/case300.m"
data = DPM.parse_file(case_path)

partition_path= "test/data/case300_3areas.csv"
DPM.assign_area!(data, partition_path)

## Settings
alpha = 1.02;
max_iteration = 1000;
tol = 1e-6;
distributed_algorithm = "ATC" # Select algorithm among APP ADMM ATC
optimizer = Gurobi.Optimizer
pf_model = "DC" # Select model among DC AC SOC SDP QC

settings = DPM.set_setting(distributed_algorithm, tol, max_iteration)

## Run distirbuted opf
subsystem = DPM.run_dopf(data, pf_model, optimizer, settings, alpha, verbose = true)

DPM.compare_solution(data, subsystem, pf_model, optimizer)
