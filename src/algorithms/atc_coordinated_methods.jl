###############################################################################
#             Build methods for ATC algorithm with coordinator                #
###############################################################################

"""
ATC algorithm module contains build and update methods
"""
module atc_coordinated_methods
using ..PowerModelsADA

"solve distributed OPF using ATC algorithm with central coordinator"
function solve_method(data, model_type::DataType, optimizer; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data=["solution", "mismatch"], print_level::Int64=1, alpha::Real=1.05, beta::Real=1.0, initialization_method::String="flat")
    solve_dopf_coordinated(data, model_type, optimizer, atc_coordinated_methods; mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, save_data=save_data, print_level=print_level, alpha=alpha, beta=beta, initialization_method=initialization_method)
end

"initialize the ATC algorithm local area"
function initialize_method_local(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    area_id = get_area_id(data)
    initialization_method = get(kwargs, :initialization_method, "flat")

    # primal and dual shared variables
    data["shared_variable"] = initialize_shared_variable(data, model_type, area_id, [0], "shared_variable", initialization_method)

    data["received_variable"] = initialize_shared_variable(data, model_type, area_id, [0], "received_variable", initialization_method)

    data["dual_variable"] = initialize_shared_variable(data, model_type, area_id ,[0], "dual_variable", initialization_method)

    # distributed algorithm settings
    initialize_dopf!(data, model_type; kwargs...)
 
    # initialize ATC parameters
    data["parameter"] = Dict( 
        "alpha" => get(kwargs, :alpha, 1.05),
        "beta" => get(kwargs, :beta, 1),
        "beta_max" => get(kwargs, :beta_max, 1e6))

end

"initialize the ATC algorithm coordinator"
function initialize_method_coordinator(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
    initialization_method = get(kwargs, :initialization_method, "flat")

    # initialize primal and dual shared variables
    data["shared_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", initialization_method)

    data["received_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "received_variable", initialization_method)

    data["dual_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "dual_variable", initialization_method)

    # distributed algorithm settings
    initialize_dopf!(data, model_type; kwargs...)

    # initialize ATC parameters
    data["parameter"] = Dict( 
        "alpha" => get(kwargs, :alpha, 1.05),
        "beta" => get(kwargs, :beta, 1),
        "beta_max" => get(kwargs, :beta_max, 1e6))

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
    beta = pm.data["parameter"]["beta"]

    ## data
    shared_variable_received = pm.data["received_variable"]
    dual_variable = pm.data["dual_variable"]

    ## objective function
    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PowerModelsADA._var(pm, variable, idx)
                v_central = shared_variable_received[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]
 
                objective += (beta * (v - v_central))^2 + v_dual * (v - v_central)
            end
        end
    end

    return objective
end

"ATC algorithm objective local area"
objective_atc_coordinator(pm::AbstractPowerModel) = objective_atc_local(pm)

"update the ATC algorithm coordinator data after each iteration"
function update_method_local(data::Dict{String, <:Any})

    ## ATC parameters
    alpha = data["parameter"]["alpha"]
    beta  = data["parameter"]["beta"]
    beta_max  = data["parameter"]["beta_max"]

    ## data
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]
    dual_variable = data["dual_variable"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_local = shared_variable_local[area][variable][idx]
                v_central =  shared_variable_received[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["dual_variable"][area][variable][idx]= v_dual + 2 * beta^2 * (v_local - v_central)
            end
        end
    end

    ## update ATC parameter
    if beta < beta_max
        data["parameter"]["beta"] *= alpha
    end

    calc_mismatch!(data)
    save_solution!(data)
    update_flag_convergence!(data)
    update_iteration!(data)
end

"update the ATC algorithm coordinator data after each iteration"
update_method_coordinator(data::Dict{String, <:Any}) = update_method_local(data)

post_processors_local = [update_solution!, update_shared_variable!]

post_processors_coordinator = [update_solution!, update_shared_variable!]

end

"""
    solve_dopf_atc_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, 
    max_iteration::Int64=1000, print_level = true, alpha::Real=1000)

Solve the distributed OPF problem using ATC algorithm with central coordinator.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- print_level::Int64=1 : print mismatch after each iteration and result summary 
- alpha::Real=1.05 : algorithm parameters
- beta::Real=1.0 : algorithm parameters
"""
solve_dopf_atc_coordinated = atc_coordinated_methods.solve_method

# export the algorithm methods module and solve method
export atc_coordinated_methods, solve_dopf_atc_coordinated