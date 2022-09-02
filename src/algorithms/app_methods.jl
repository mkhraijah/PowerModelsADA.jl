###############################################################################
#                     Build methods for APP algorithm                        #
###############################################################################

"""
    solve_dopf_app(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, alpha::Real=1000, beta::Real, gamma::Real)

Solve the distributed OPF problem using APP algorithm.
"""
function solve_dopf_app(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, alpha::Real=1000, beta::Real=2alpha, gamma::Real=alpha)
    solve_dopf(data, model_type, optimizer, initialize_dopf_app!, update_app!, build_dopf_app ; tol=tol , max_iteration=max_iteration, verbose = verbose, alpha=alpha, beta=beta, gamma=gamma)
end

"inilitlize the APP algorithm"
function initialize_dopf_app!(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)

    # initiate APP parameters
    alpha = get(kwargs, :alpha, 1000)
    beta = get(kwargs, :beta, 2*alpha)
    gamma = get(kwargs, :gamma, alpha)

    data["alpha"] = alpha
    data["beta"] = beta
    data["gamma"] = gamma

end

"build PowerModel using APP algorithm"
function build_dopf_app(pm::AbstractPowerModel)

    # define variables
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    # define constraints
    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
    end

    for i in ids(pm, :dcline)
        _PM.constraint_dcline_power_losses(pm, i)
    end

    objective_min_fuel_and_consensus!(pm, objective_app!)
end

"set the APP algorithm objective"
function objective_app!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.data["alpha"]
    beta  = pm.data["beta"]
    gamma = pm.data["gamma"]

    ## data
    area_id = string(get_area_id(pm))
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## objective function

    objective = JuMP.objective_function(pm.model)
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = _var(pm, variable, idx)
                v_neighbor = primal_variable[area][variable][idx]
                v_local = primal_variable[area_id][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                objective += beta/2 *  (v - v_local)^2 + gamma * v * (v_local - v_neighbor) + v * v_dual
            end
        end
    end

    JuMP.@objective(pm.model, Min,  objective)
end

"update the APP algorithm before each iteration"
function update_app!(data::Dict{String, <:Any})

    ## APP parameters
    alpha = data["alpha"]

    ## data
    area_id = string(get_area_id(data))
    primal_variable = data["shared_primal"]
    dual_variable = data["shared_dual"]

    ## update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_neighbor = primal_variable[area][variable][idx]
                v_local = primal_variable[area_id][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["shared_dual"][area][variable][idx] = v_dual  + alpha * (v_local - v_neighbor)
            end
        end
    end

end
