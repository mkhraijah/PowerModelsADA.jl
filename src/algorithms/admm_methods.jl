###############################################################################
#                     Build methods for ADMM algorithm                        #
###############################################################################

"""
    solve_dopf_admm(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose = true, print_optimizer_info::Bool=false, alpha::Real=1000)

Solve the distributed OPF problem using ADMM algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Bool=true : print mismatch after each iteration and result summary 
- print_optimizer_info::Bool=false : print local optimization info from the solver
- alpha = algorithm parameters
"""
function solve_dopf_admm(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Bool=true, print_optimizer_info::Bool=false, alpha::Real=1000)

    solve_dopf(data, model_type, optimizer, admm_methods; 
    mismatch_method=mismatch_method, tol=tol , max_iteration=max_iteration, 
    verbose=verbose, print_optimizer_info=print_optimizer_info, alpha=alpha)
end

export solve_dopf_admm

"""
ADMM algorithm module contians build and update methods
"""
module admm_methods
using ..PMADA

"inilitlize the ADMM algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::DataType; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)

    # initiate ADMM parameters
    alpha = get(kwargs, :alpha, 1000)
    data["alpha"] = alpha
end

"build PowerModel using ADMM algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_admm)
end

"ADMM algorithm objective"
function objective_admm(pm::AbstractPowerModel)

    ## ADMM parameters
    alpha = pm.data["alpha"]

    ## data
    area_id = string(get_area_id(pm))
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## objective function
    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PMADA._var(pm, variable, idx)
                v_central = (primal_variable[area_id][variable][idx] + primal_variable[area][variable][idx])/2
                v_dual = dual_variable[area][variable][idx]
 
                objective += alpha/2 * (v - v_central)^2 + v_dual * (v - v_central)
            end
        end
    end

    return objective
end


"update the ADMM algorithm before each iteration"
function update_method(data::Dict{String, <:Any})

    ## ADMM parameters
    alpha = data["alpha"]

    ## data
    area_id = string(get_area_id(data))
    primal_variable = data["shared_primal"]
    dual_variable = data["shared_dual"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_primal = primal_variable[area_id][variable][idx]
                v_central = (primal_variable[area_id][variable][idx] + primal_variable[area][variable][idx])/2
                v_dual = dual_variable[area][variable][idx]

                data["shared_dual"][area][variable][idx]= v_dual  + alpha * (v_primal - v_central)
            end
        end
    end
end
end