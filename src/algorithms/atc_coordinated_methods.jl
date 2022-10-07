###############################################################################
#             Build methods for ATC algorithm with coordinator                #
###############################################################################

"""
ATC algorithm module contians build and update methods
"""
module atc_coordinated_methods
using ..PMADA

"solve distributed OPF using ATC algorithm with central coordinator"
function solve_method(data, model_type::DataType, 
    optimizer; mismatch_method::String="norm", tol::Float64=1e-4, 
    max_iteration::Int64=1000, save_data=["solution", "mismatch"], 
    verbose::Int64=1, alpha::Real=1.05, beta::Real=1.0)
    
    solve_dopf_coordinated(data, model_type, optimizer, atc_coordinated_methods; 
    mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, 
    save_data=save_data, verbose=verbose, alpha=alpha, beta=beta)
end

"inilitlize the ATC algorithm local area"
function initialize_method_local(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    # primal and dual shared variables
    initialize_shared_variable_local!(data, model_type)

    initialize_dual_variable_local!(data, model_type)

    # solution dictionary
    initialize_solution!(data, model_type)

    # distributed algorithm parameters
    initialize_dopf_parameters!(data; kwargs...)
 
    # atc parameters
    data["alpha"] = alpha = get(kwargs, :alpha, 1.05)
    data["beta"] = beta = get(kwargs, :beta, 1)
end

"inilitlize the ATC algorithm coordinator"
function initialize_method_coordinator(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    # initiate primal and dual shared variables
    initialize_shared_variable_coordinator!(data, model_type)

    initialize_dual_variable_coordinator!(data, model_type)

    # initiate solution dictionary
    initialize_solution!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; kwargs...)

    # initiate atc parameters
    data["alpha"] = get(kwargs, :alpha, 1.05)
    data["beta"] = get(kwargs, :beta, 1)
end

"build PowerModel for ATC algorithm local area"
function build_method_local(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_atc_local)
end

"build PowerModel for ATC algorithm coordinator"
function build_method_coordinator(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_atc_coordinator)
end

"ATC algorithm objective coordinator"
function objective_atc_local(pm::AbstractPowerModel)

    ## atc parameters
    beta = pm.data["beta"]

    ## data
    shared_variable = pm.data["shared_variable"]
    dual_variable = pm.data["dual_variable"]

    ## objective function
    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PMADA._var(pm, variable, idx)
                v_central = shared_variable[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]
 
                objective += (beta * (v - v_central))^2 + v_dual * (v - v_central)
            end
        end
    end

    return objective
end

"ATC algorithm objective local area"
objective_atc_coordinator(pm::AbstractPowerModel) = objective_atc_local(pm)

"update the ATC algorithm coordinator data before each iteration"
function update_method_local_before(data::Dict{String, <:Any})

    ## ATC parameters
    alpha = data["alpha"]
    beta  = data["beta"]

    ## data
    area_id = string(get_area_id(data))
    shared_variable = data["shared_variable"]
    dual_variable = data["dual_variable"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_primal = shared_variable[area_id][variable][idx]
                v_local =  shared_variable[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["dual_variable"][area][variable][idx]= v_dual  + 2 * beta^2 * (v_primal - v_local)
            end
        end
    end
    ## update ATC parameter
    data["beta"] *= alpha
end

"update the ATC algorithm local area data before each iteration"
update_method_coordinator_before(data::Dict{String, <:Any}) = update_method_local_before(data)

"update the ATC algorithm local area data between each iteration"
function update_method_local_between(data::Dict{String, <:Any})

end

"update the ATC algorithm coordinator data between each iteration"
update_method_coordinator_between(data::Dict{String, <:Any}) = update_method_local_between(data)

"update the ATC algorithm coordinator data after each iteration"
function update_method_local_after(data::Dict{String, <:Any})
    
    save_solution!(data)
    calc_mismatch!(data)
    update_flag_convergance!(data)
    update_iteration!(data)
end

"update the ATC algorithm coordinator data after each iteration"
update_method_coordinator_after(data::Dict{String, <:Any}) = update_method_local_after(data)
end

"""
    solve_dopf_atc_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, 
    max_iteration::Int64=1000, verbose = true, alpha::Real=1000)

Solve the distributed OPF problem using ATC algorithm with central coordinator.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Int64=1 : print mismatch after each iteration and result summary 
- alpha::Real=1.05 : algorithm parameters
- beta::Real=1.0 : algorithm parameters
"""
solve_dopf_atc_coordinated = atc_coordinated_methods.solve_method

# export the algorithm methods module and call method
export atc_coordinated_methods, solve_dopf_atc_coordinated