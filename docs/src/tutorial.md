# Tutorial 

PMADA solve the OPF problem using either pre-defined distributed algorithm or user-defined algorithm. This page shows example of solving the OPF problem using the pre-defined algorithms and how to define a new alternating distributed algorithm. 

## Run Distributed Algorithm
To solve the OPF algorithm, we need first to import `PMADA` package and an optimization solver. In this case we use `Ipopt` a NLP solver. You can install the solver using `using Pkg, Pkg.add("Ipopt")`. Then run the following code while you are inside the PMADA package directory. 
```julia
## Import package
using PMADA
using Ipopt 
```

Next, we need to upload a test case. We will use IEEE 14-bus system in `/test/data/` folder in matpower format. The file can be loaded using `parse_file` from PowerModels package. The test system needs to be divided into different areas. This can be check by looking into `data["bus"][bus_id]["area"]`. If not, you can use `partition_system!` function to divide the system or load a csv file that constains the buses and area of each bus.  

```julia

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
parition_file_path = "test/data/case14_2areas.csv"
data = parse_file(case_path)
assign_area!(data, parition_file_path)

```
Now, the case study is loaded and ready to solve the OPF problem using distirbuted algorithms. We first need to define parameters, load solver, and select power flow formulaiton `model_type` as follow: 


```julia

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
alpha = 1000
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

##  Power Flow Model selection
model_type = DCPPowerModel

```

To solve the OPF problem using APP algorithm using the call function, we use the following code: 

```julia
##  Distributed algorithm
data_area = solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = false, alpha=alpha);
error_admm = compare_solution(data, data_area, model_type, optimizer)
```

To solve the OPF problem using APP algorithm, we use the following code: 

```julia
##  Distributed algorithm
data_area = solve_dopf_app(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, verbose = false, alpha=alpha);
error_admm = compare_solution(data, data_area, model_type, optimizer)
```

Another way to solve the OPF problem in more controlled method. We can use the following code: 

```julia

## Import package
using PMADA
using Ipopt 

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
alpha = 1000
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
parition_file_path = "test/data/case14_2areas.csv"
data = parse_file(case_path)
assign_area!(data, parition_file_path)

## obtain areas idx
areas_id = get_areas_id(data)

## decompose the system into subsystems
data_area = decompose_system(data)

## initilize distributed PowerModels parameters
for i in areas_id
    admm_methods.initialize_method(data_area[i], model_type; tol=tol, max_iteration=max_iteration, kwargs...)
end

## initialaize the algorithms global counters
iteration = 0
flag_convergance = false

## start iteration
while iteration < max_iteration && flag_convergance == false

    
    ## solve local problem and update solution
    for i in areas_id
        admm_methods.update_method(data_area[i])
        solve_local!(data_area[i], model_type, optimizer, admm_methods.build_method)
    end

    ## share solution with neighbors, the shared data is first obtained to facilitate distributed implementation  
    for i in areas_id # sender subsystem
        for j in areas_id # receiver subsystem
            if i != j && string(i) in keys(data_area[j]["shared_primal"])
                shared_data = send_shared_data(i, j, data_area[i])
                receive_shared_data!(i, shared_data, data_area[j])
            end
        end
    end

    ## calculate mismatches and update convergance flags
    for i in areas_id
        calc_mismatch!(data_area[i])
        update_flag_convergance!(data_area[i], tol)
        update_iteration!(data_area[i])
    end

    ## check global convergance and update iteration counters
    flag_convergance = update_global_flag_convergance(data_area)
    iteration += 1
end

```