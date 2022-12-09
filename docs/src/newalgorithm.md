# User-Defined Algorithm

```@meta
CurrentModule = PowerModelsADA
```

To define a new algorithm, we need to define a module for the new algorithm that contains the main solve function in addition to three algorithm-specific functions. The three algorithm-specific are: initialize, build, and update. You can follow the exmaple in the [template file](https://github.com/mkhraijah/PowerModelsADA.jl/blob/main/example/template.jl).

The module of `xx` algorithm should be defined and exported as `xx_methods` as follows:

```julia
"""
template for xx distributed algorithm
"""
module xx_methods
using ..PowerModelsADA

### functions ###

# export the algorithm methods module and call method
export xx_methods, solve_dopf_xx
```

The solve function is the main method to use the `xx` algorithm. The function takes the data, power flow formulation (`model_type`), JuMP solver object, and algorithm's parameters as required. The solve function should use the pre-defined algorithm flow as follows:  

```julia
"solve distributed OPF using XX algorithm"
function solve_method(data, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Int64=1, parameters...)

    solve_dopf(data, model_type, optimizer, xx_methods; 
    mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, 
    verbose=verbose, parameters...)
end
```

The first algorithm-specific function is the initialize function. The function takes the area data file and adds to it the required parameters, counters, and shared variables. There are multiple built-in functions in `PowerModelsADA` that can be used to define the shared and received variables, as well as the dual variables. Note that the initialization function should include the `initialize_dopf!` to define the counters and convergence flags.

```julia
"initialize the XX algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    data["shared_variable"] = Dict(to_area=> variable_name=>value)
    data["received_variable"] = Dict(from_area=> variable_name=>value)

    # distributed algorithm settings
    initialize_dopf!(data, model_type; kwargs...)

    # xx parameters
    data["parameter"] = Dict("alpha"=> get(kwargs, :alpha, 1000))

end
```

The second function is the build function, which builds the `PowerModels` object of the subproblem. The subproblems typically have the same variables and constraints as the central OPF problem and differ in the objective functions. To build a subproblem with the same variables and constraints as the central OPF problem with a specific objective function, we need to define the objective function using the template shown below. The objective function definition takes the `PowerModels` object and returns a `JuMP` expression. You can use the internal helper function `_var` to obtain the `JuMP` model variables' object defined in the `PowerModels` object.

```julia
"build PowerModel using xx algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_function)
end

"set the xx algorithm objective"
function objective_function(pm::AbstractPowerModel)

    # to get the JuMP object of the active power of generator 1 use:
    pg1 = _var(pm, :pg, 1)

    ###
    objective = pg1
    ###
    return objective
end
```

```@docs
_var
```

The last function is to update the area dictionary after communicating the shared variables results with other areas. 

```julia
"update the xx algorithm before each iteration"
function update_method(data::Dict{String, <:Any})

    ###

    ###
end

```

This is a general way to define a distributed algorithm that is fully distributed with the same main algorithm flow as the pre-defined algorithms. For other algorithm flows, the solve function needs to be defined fully instead of using the pre-define function `solve_dopf`.
