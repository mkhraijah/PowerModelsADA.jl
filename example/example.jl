## Import package
using PowerModelsADA
using Ipopt

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
parition_file_path = "test/data/case14_2areas.csv"
data = parse_file(case_path)
assign_area!(data, parition_file_path)

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-2
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
model_type = ACPPowerModel

## Distributed algorithm
## ADMM with fully distributed structure
data_area = solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, print_level = 1, alpha=1000, save_data=["solution", "mismatch"])
error_admm = compare_solution(data, data_area, model_type, optimizer)

## APP with fully distributed structure
data_area = solve_dopf_app(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, print_level = 1, alpha=1000, save_data=["solution", "mismatch"])
error_app = compare_solution(data, data_area, model_type, optimizer)

## ATC with fully distributed structure
data_area = solve_dopf_atc(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, print_level = 1, alpha=1.1)
error_atc = compare_solution(data, data_area, model_type, optimizer)

## ADMM with central coordinator structure
data_area = solve_dopf_admm_coordinated(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, print_level = 1, alpha = 1000);
error_admm = compare_solution(data, data_area, model_type, optimizer)

## ATC with central coordinator structure
data_area = solve_dopf_atc_coordinated(data, model_type, optimizer; max_iteration=max_iteration, print_level = 1, alpha = 1.05)
error_atc = compare_solution(data, data_area, model_type, optimizer)

## ALADIN with central coordinator structure
sigma = Dict{String, Real}("va" => 10, "vm" => 5, "pf" => 1, "pt" => 1, "qf" => 1, "qt" => 1, "pg" => 1, "qg" => 1)
data_area = solve_dopf_aladin_coordinated(data, model_type, optimizer; tol=tol, max_iteration=max_iteration, print_level=1, p=100, mu=1000, r_p=1.5, r_mu=2, q_gamma=0, sigma=sigma)
error_aladin = compare_solution(data, data_area, model_type, optimizer)