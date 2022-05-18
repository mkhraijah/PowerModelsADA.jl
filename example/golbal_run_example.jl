
path = joinpath(@__DIR__,"..")
cd(path)
using Pkg
Pkg.activate(".")
using Ipopt
using JuMP

## Call Packages
include("../src/DistributedPowerModels.jl")
DPM = DistributedPowerModels

## Read case function with partition file return dictionary of the case
case_path = "test/data/case14.m"
data = DPM.parse_file(case_path)

partition_path= "test/data/case300_3areas.csv"
DPM.assign_area!(data, partition_path)

## Settings
alpha = 100000;
max_iteration = 1000;
tol = 1e-2;
optimizer = Juniper.Optimizer("nl_solver"=>Ipopt.Optimizer)

nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
optimizer= optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver)

optimizer= SCS.Optimizer
pf_model = "SDP" # Select model among DC AC SOC SDP QC
build_method = DPM.build_dopf_admm


## Run distirbuted opf
data_area = DPM.run_dopf(data_WB5, pf_model, build_method, optimizer, alpha, verbose = true)

DPM.compare_solution(data, data_area, pf_model, optimizer)
