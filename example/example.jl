## Import package
include("../src/PMADA.jl")
using .PMADA 
using Ipopt


## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
parition_file_path = "test/data/case14_2areas.csv"
data = PMADA.parse_file(case_path)
PMADA.assign_area!(data, parition_file_path)

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
optimizer = PMADA.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

##  Power Flow Model selection
model_type = PMADA.DCPPowerModel

##  Distributed algorithm
data_area = PMADA.solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = 1, alpha=10000);
error_admm = PMADA.compare_solution(data, data_area, model_type, optimizer)
data_area = solve_dopf_app(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, verbose = 1, alpha=10000);
error_app = compare_solution(data, data_area, model_type, optimizer)
data_area = solve_dopf_atc(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = 1, alpha=1.05);
error_atc = compare_solution(data, data_area, model_type, optimizer)


data_area = solve_dopf_admm_coordinated(data, model_type, optimizer; max_iteration=max_iteration, verbose = 1, alpha = 100000);
error_admm = compare_solution(data, data_area, model_type, optimizer)
data_area = solve_dopf_atc_coordinated(data, model_type, optimizer; max_iteration=max_iteration, verbose = 1, alpha = 1.05);
error_atc = compare_solution(data, data_area, model_type, optimizer)