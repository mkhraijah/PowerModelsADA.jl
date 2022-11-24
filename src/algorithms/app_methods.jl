###############################################################################
#                     Build methods for APP algorithm                        #
###############################################################################

"""
APP algorithm module contains build and update methods
"""
module app_methods
using ..PowerModelsADA

"solve distributed OPF using APP algorithm"
function solve_method(data, model_type::DataType, optimizer; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data=["solution"], print_level::Int64=1, alpha::Real=1000, beta::Real=2alpha, gamma::Real=alpha, initialization_method::String="flat")

    solve_dopf(data, model_type, optimizer, app_methods; mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, save_data=save_data, print_level=print_level, alpha=alpha, beta=beta, gamma=gamma, initialization_method=initialization_method)
end

"initialize the APP algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::Type; kwargs...)

    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id
    initialization_method = get(kwargs, :initialization_method, "flat")

    # primal and dual shared variables
    data["shared_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", initialization_method)

    data["received_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "received_variable", initialization_method)

    data["dual_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "dual_variable", initialization_method)

    # distributed algorithm settings
    initialize_dopf!(data, model_type; kwargs...)

    # initialize APP parameters
    data["parameter"] = Dict( 
        "alpha" => get(kwargs, :alpha, 1000),
        "beta" => get(kwargs, :beta, 2*get(kwargs, :alpha, 1000)),
        "gamma" => get(kwargs, :gamma, get(kwargs, :alpha, 1000)))

end

"build PowerModel using APP algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_app)
end

"set the APP algorithm objective"
function objective_app(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.data["parameter"]["alpha"]
    beta  = pm.data["parameter"]["beta"]
    gamma = pm.data["parameter"]["gamma"]

    ## data
    shared_variable_local = pm.data["shared_variable"]
    shared_variable_received = pm.data["received_variable"]
    dual_variable = pm.data["dual_variable"]

    ## objective function

    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PowerModelsADA._var(pm, variable, idx)
                v_neighbor = shared_variable_received[area][variable][idx]
                v_local = shared_variable_local[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                objective += beta/2 *  (v - v_local)^2 + gamma * v * (v_local - v_neighbor) + v * v_dual
            end
        end
    end

    return objective
end

"update the APP algorithm data after each iteration"
function update_method(data::Dict{String, <:Any})
    
    ## APP parameters
    alpha = data["parameter"]["alpha"]

    ## data
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]
    dual_variable = data["dual_variable"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_neighbor = shared_variable_received[area][variable][idx]
                v_local = shared_variable_local[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["dual_variable"][area][variable][idx] = v_dual  + alpha * (v_local - v_neighbor)
            end
        end
    end

    calc_mismatch!(data)
    update_flag_convergence!(data)
    save_solution!(data)
    update_iteration!(data)
end

post_processors = [update_solution!, update_shared_variable!]
end

"""
    solve_dopf_app(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm",tol::Float64=1e-4, max_iteration::Int64=1000, 
    print_level::Int64=1, alpha::Real=1000, beta::Real, gamma::Real)

Solve the distributed OPF problem using APP algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- print_level::Int64=1 : print mismatch after each iteration and result summary 
- alpha::Real= 1000 : algorithm parameters
- beta::Real= 2alpha : algorithm parameters
- gamma::Real= alpha : algorithm parameters
"""
solve_dopf_app = app_methods.solve_method

# export the algorithm methods module and solve method
export app_methods, solve_dopf_app