# Tutorial

`PowerModelsADA` solves OPF problems using either a pre-defined distributed algorithm or a user-defined algorithm. This page shows examples of solving the OPF problem using the pre-defined algorithms and how to define a new alternating distributed algorithm.

The distributed algorithm-specific functions are stored in modules. Each module contains at least three main functions: initialize, build, and update functions. Each module also contains a solve function that solves the OPF by passing the case, solver, and power flow model.

The distributed algorithm module and solve function are:

- ADMM: modules: `admm_methods` and `admm_coordinated_methods`. solve functions: `solve_dopf_admm` and `solve_dopf_admm_coordinated`
- ATC: modules: `atc_methods` and `atc_coordinated_methods`. solve functions: `solve_dopf_atc` and `solve_dopf_atc_coordinated`
- APP: modules: `app_methods`. solve functions: `solve_dopf_app`
- ALADIN: modules: `aladin_coordinated_methods`. solve function: `solve_dopf_aladin_coordinated`
- Adaptive ADMM: modules: `adaptive_admm_methods` and `adaptive_admm_coordinated_methods`. solve functions: `solve_dopf_adaptive_admm` and `solve_dopf_adaptive_admm_coordinated`

## Run Distributed Algorithm

To solve the OPF problem, we need first to import the `PowerModelsADA` package and an optimization solver. In this case we use the NLP solver `Ipopt`. You can install the solver using `using Pkg; Pkg.add("Ipopt")`. Then run the following code:

```julia
## Import package
using PowerModelsADA
using Ipopt 
```

Next, we need to upload a test case. We will use IEEE 14-bus system in `/test/data/` folder in MATPOWER format. The file can be loaded using `parse_file` from `PowerModels` package. The test system needs to be divided into multiple distinct areas. This can be checked by looking into `data["bus"][bus_id]["area"]`.

```julia
## Read case with partition file and return dictionary of the partitioned case
case_path = "test/data/case14.m"
partition_file_path = "test/data/case14_2areas.csv"
data = parse_file(case_path)
assign_area!(data, partition_file_path)
```

Now, the case study is loaded and ready to be used to solve the OPF problem using distributed algorithms. We first need to define parameters, load the solver, and select a power flow formulation `model_type` as follows:

```julia

## Settings and optimizer initiation
max_iteration = 1000
tol = 1e-4
alpha = 1000
optimizer = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
```

PowerModelsADA supports the following power flow models:

### Exact power flow

```julia
model_type = ACPPowerModel # AC power flow model with polar bus voltage variables.
model_type = ACRPowerModel # AC power flow model with rectangular bus voltage variables.
```

### Approximations

```julia
model_type = DCPPowerModel # Linearized 'DC' power flow model.
model_type = LPACCPowerModel # LP AC power flow approximation.
```

### Convex relaxations

```julia
model_type = SOCWRPowerModel # Second-order cone relaxation of bus injection model of AC power flow.
model_type = QCRMPowerModel # Quadratic-Convex relaxation of the AC power flow.
model_type = SDPWRMPowerModel # Semidefinite relaxation of AC power flow.
model_type = SparseSDPWRMPowerModel # Sparsity-exploiting semidefinite relaxation of AC power flow.
```

To solve the OPF problem using ADMM algorithm using the solve function, we use the following:

```julia
data_area = solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, alpha=alpha)

```

To use multiprocessing features, we need to use the Distributed library, add processors, and upload the PowerModelsADA and the solver packages to the processors. For the best performance, the number of processors should be equal to the number of areas. The code becomes as follows:

```julia
using Distributed 
num_area = 4 # change the number to be equal to number of areas
addprocs(num_area, exeflags="--project")
@everywhere using PowerModelsADA
@everywhere using Ipopt
data_area = solve_dopf_admm(data, model_type, optimizer, tol=tol, max_iteration=max_iteration, alpha=alpha, multiprocessors=true)
```

To compare the distributed algorithm objective function value with the central OPF, use `compare_solution` to get the absolute value of the relative error.

```julia
optimality_gap = compare_solution(data, data_area, model_type, optimizer)
```

PowerModelsADA also provides the flexibility for more granular control of the distributed algorithm. We can use the following code to initialize the distributed algorithm (we use ADMM in this example).

```julia

## define parameters and power flow model
max_iteration = 1000
tol = 1e-4
alpha = 1000
model_type = DCPPowerModel

## obtain areas idx
areas_id = get_areas_id(data)

## decompose the system into subsystems
data_area = decompose_system(data)

## initialize parameters using the algorithm-specific initialize function
for i in areas_id
    admm_methods.initialize_method(data_area[i], model_type; tol=tol, max_iteration=max_iteration, alpha = alpha)
end

```

We then start the iterative process of the distributed algorithm using while loop with a pre-define termination criteria as follows:

```julia
## initialize global counters
iteration = 0
flag_convergence = false

## start iteration
while iteration < max_iteration && flag_convergence == false

    ## solve local problem and update solution
    for i in areas_id
        result = solve_model(data_area[i], model_type, optimizer, admm_methods.build_method, solution_processors=admm_methods.post_processors)
        update_data!(data_area[i], result["solution"])
    end

    ## share solution with neighbors
    for i in areas_id # sender subsystem
        for j in data_area[i]["neighbors"] # receiver subsystem
            shared_data = prepare_shared_data(data_area[i], j)
            receive_shared_data!(data_area[j], deepcopy(shared_data), i)
        end
    end

    # calculate mismatches and update convergence flags
    for i in areas_id
        dopf_method.update_method(data_area[i])
    end

    ## check global convergence and update iteration counters
    flag_convergence = update_global_flag_convergence(data_area)
    iteration += 1
end

```
