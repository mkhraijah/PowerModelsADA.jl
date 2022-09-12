###############################################################################
#                     Build methods for ATC algorithm                        #
###############################################################################

"""
    solve_dopf_atc(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose = true, print_optimizer_info::Bool=false, alpha::Real=1000, beta::Real = 1)

Solve the distributed OPF problem using ATC algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Bool=true : print mismatch after each iteration and result summary 
- print_optimizer_info::Bool=false : print local optimization info from the solver
- alpha::Real=1.05 : algorithm parameters
- beta::Real=1.0 : algorithm parameters
"""
function solve_dopf_atc(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Bool=true, print_optimizer_info::Bool=false, alpha::Real=1000, beta::Real=1)

    solve_dopf(data, model_type, optimizer, atc_methods; 
    mismatch_method=mismatch_method, tol=tol , max_iteration=max_iteration, 
    verbose=verbose, print_optimizer_info=print_optimizer_info, alpha=alpha, beta=beta)
end

function solve_dopf_atc(data::String, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    verbose::Bool=true, print_optimizer_info::Bool=false, alpha::Real=1000, beta::Real=1)

    solve_dopf(data, model_type, optimizer, atc_methods; 
    mismatch_method=mismatch_method, tol=tol , max_iteration=max_iteration, 
    verbose=verbose, print_optimizer_info=print_optimizer_info, alpha=alpha, beta=beta)
end

export solve_dopf_atc

"""
ATC algorithm module contians build and update methods
"""
module atc_methods
using ..PMADA

"redefine the main call method inside the module"
solve_method = solve_dopf_atc

"inilitlize the ATC algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::DataType; 
    tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)

    # initiate ATC parameters
    alpha = get(kwargs, :alpha, 1000)
    beta = get(kwargs, :beta, 1)

    data["alpha"] = alpha
    data["beta"] = beta
   
end

"build method for Distributed PowerModel using ATC algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_atc)
end

"set the ATC algorithm objective"
function objective_atc(pm::AbstractPowerModel)

    ## ATC parameters
    beta = pm.data["beta"]

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

                objective += (beta * (v - v_central))^2 + v_dual * (v - v_central)
            end
        end
    end

    return objective
end


"update the ATC algorithm before each iteration"
function update_method(data::Dict{String, <:Any})

    ## ATC parameters
    alpha = data["alpha"]
    beta  = data["beta"]

    ## data
    area_id = string(get_area_id(data))
    primal_variable = data["shared_primal"]
    dual_variable = data["shared_dual"]

    ## update dual variable
    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_primal = primal_variable[area_id][variable][idx]
                v_central = (primal_variable[area_id][variable][idx] + primal_variable[area][variable][idx])/2
                v_dual = dual_variable[area][variable][idx]

                data["shared_dual"][area][variable][idx] = v_dual  + 2 * beta^2 * (v_primal - v_central)
            end
        end
    end

    ## update ATC parameter
    data["beta"] *= alpha

end
end

# export the algorithm methods module and call method
export atc_methods