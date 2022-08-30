## Import package
using PMADA
using Ipopt
using JuMP



## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case300.m"
data = PMADA.parse_file(case_path)


## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
# optimizer = JuMP.optimizer_with_attributes(Hypatia.Optimizer, "verbose"=>0)

##  Power Flow Model selection

model_type = PMADA.SparseSDPWRMPowerModel

##  Distributed algorithm settings

data_area = PMADA.solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=1000, verbose = true, alpha=1000)
error = PMADA.compare_solution(data, data_area, model_type, optimizer)
data_area = PMADA.solve_dopf_app(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, verbose = true, alpha=1000)
error = PMADA.compare_solution(data, data_area, model_type, optimizer)
data_area = PMADA.solve_dopf_atc(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = true, alpha=1.05)
error = PMADA.compare_solution(data, data_area, model_type, optimizer)