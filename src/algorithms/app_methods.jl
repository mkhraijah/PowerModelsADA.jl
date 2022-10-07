###############################################################################
#                     Build methods for APP algorithm                        #
###############################################################################

"""
APP algorithm module contians build and update methods
"""
module app_methods
using ..PMADA

"solve distributed OPF using APP algorithm"
function solve_method(data, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    save_data=["solution"], verbose::Int64=1, 
    alpha::Real=1000, beta::Real=2alpha, gamma::Real=alpha)

    solve_dopf(data, model_type, optimizer, app_methods; 
    mismatch_method=mismatch_method, tol=tol, max_iteration=max_iteration, 
    save_data=save_data, verbose=verbose, alpha=alpha, beta=beta, gamma=gamma)
end

"inilitlize the APP algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::Type; kwargs...)

    # primal and dual shared variables
    initialize_shared_variable!(data, model_type)

    initialize_dual_variable!(data, model_type)

    # solution dictionary
    initialize_solution!(data, model_type)

    # distributed algorithm parameters
    initialize_dopf_parameters!(data; kwargs...)

    # APP parameters
    data["alpha"] = get(kwargs, :alpha, 1000)
    data["beta"] = get(kwargs, :beta, 2*data["alpha"])
    data["gamma"] = get(kwargs, :gamma, data["alpha"])

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
    alpha = pm.data["alpha"]
    beta  = pm.data["beta"]
    gamma = pm.data["gamma"]

    ## data
    area_id = string(get_area_id(pm))
    shared_variable = pm.data["shared_variable"]
    dual_variable = pm.data["dual_variable"]

    ## objective function

    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PMADA._var(pm, variable, idx)
                v_neighbor = shared_variable[area][variable][idx]
                v_local = shared_variable[area_id][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                objective += beta/2 *  (v - v_local)^2 + gamma * v * (v_local - v_neighbor) + v * v_dual
            end
        end
    end

    return objective
end

"update the APP algorithm before each iteration"
function update_method_before(data::Dict{String, <:Any})

    ## APP parameters
    alpha = data["alpha"]

    ## data
    area_id = string(get_area_id(data))
    shared_variable = data["shared_variable"]
    dual_variable = data["dual_variable"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_neighbor = shared_variable[area][variable][idx]
                v_local = shared_variable[area_id][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["dual_variable"][area][variable][idx] = v_dual  + alpha * (v_local - v_neighbor)
            end
        end
    end
end

"update the APP algorithm data after each iteration"
function update_method_after(data::Dict{String, <:Any})
    
    save_solution!(data)
    calc_mismatch!(data)
    update_flag_convergance!(data)
    update_iteration!(data)
end


end

"""
    solve_dopf_app(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm",tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Int64=1, alpha::Real=1000, beta::Real, gamma::Real)

Solve the distributed OPF problem using APP algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Int64=1 : print mismatch after each iteration and result summary 
- alpha::Real= 1000 : algorithm parameters
- beta::Real= 2alpha : algorithm parameters
- gamma::Real= alpha : algorithm parameters
"""
solve_dopf_app = app_methods.solve_method

# export the algorithm methods module and call method
export app_methods, solve_dopf_app