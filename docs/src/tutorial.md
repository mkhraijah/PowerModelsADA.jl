# Tutorial 

## Run Distributed Algorithm

'''julia
## Import package
using PMADA
using Ipopt 
'''



'''julia

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
parition_file_path = "test/data/case14_2areas.csv"
data = parse_file(case_path)
assign_area!(data, parition_file_path)

'''



'''julia

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

##  Power Flow Model selection
model_type = DCPPowerModel

'''



'''julia
##  Distributed algorithm
data_area = solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = false, alpha=1000);
error_admm = compare_solution(data, data_area, model_type, optimizer)
'''







## Run User Defined Distributed Algorithm 