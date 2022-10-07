###############################################################################
#                    Build methods for the XX algorithm                       #
###############################################################################

"""
Templet for xx distributed algorithm
"""
module xx_methods
using ..PMADA

"solve distributed OPF using XX algorithm"
function solve_method(data, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Int64=1, parameters...)

    solve_dopf(data, model_type, optimizer, xx_methods; 
    mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, 
    verbose=verbose, parameters...)
end

"inilitlize the XX algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)

    # initiate APP parameters
    parameters = get(kwargs, :parameters, 1000)
    data["parameters"] = parameters

end

"build PowerModel using xx algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_function)
end

"set the APP algorithm objective"
function objective_function(pm::AbstractPowerModel)

    ###
    objective = 0
    ###
    return objective
end

"update the xx algorithm before each iteration"
function update_method_before(data::Dict{String, <:Any})

    ###

    ###
end

"update the xx algorithm data after each iteration"
function update_method_after(data::Dict{String, <:Any})
    
    ###

    ###

    calc_mismatch!(data, data["mismatch_method"])
    update_flag_convergance!(data, data["tol"])
    update_iteration!(data)
end


end

"""
    solve_dopf_xx(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm",tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Int64=1, parameters)

Solve the distributed OPF problem using APP algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Int64=1 : print mismatch after each iteration and result summary 
"""
solve_dopf_xx = xx_methods.solve_method

# export the algorithm methods module and call method
export xx_methods, solve_dopf_xx