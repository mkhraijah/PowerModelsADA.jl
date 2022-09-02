###############################################################################
#            Build methods for ADMM algorithm with coordinator                #
###############################################################################

"""
    solve_dopf_admm(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, alpha::Real=1000)

Solve the distributed OPF problem using ADMM algorithm.
"""
function solve_dopf_admm_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Bool=true, alpha::Real=1000)
    
    data_coordinator, data_area = PMADA.solve_dopf(data, model_type, optimizer, 
    PMADA.initialize_dopf_admm_coordinator!, PMADA.initialize_dopf_admm_local!,
    PMADA.build_dopf_admm_coordinator, PMADA.build_dopf_admm_local, 
    PMADA.update_admm_coordinated!;tol=tol , max_iteration=max_iteration, verbose=verbose, alpha=alpha)
    
end


function initialize_dopf_admm_local!(data::Dict{String, <:Any}, model_type::DataType; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared_local!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)
 
    # initiate ADMM parameters
    alpha = get(kwargs, :alpha, 1000)
    data["alpha"] = alpha

end

function initialize_dopf_admm_coordinator!(data::Dict{String, <:Any}, model_type::DataType; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    initialize_variable_shared_coordinator!(data, model_type)

    # initiate distributed algorithm parameters
    initialize_dopf_parameters!(data; tol=tol, max_iteration=max_iteration)

    # initiate ADMM parameters
    alpha = get(kwargs, :alpha, 1000)
    data["alpha"] = alpha

end

"build PowerModel using ADMM algorithm"
function build_dopf_admm_local(pm::AbstractPowerModel)

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

    objective_min_fuel_and_consensus!(pm, objective_admm_coordinated!)
end

function build_dopf_admm_coordinator(pm::AbstractPowerModel)

    # define variables
    _PM.variable_bus_voltage(pm)
    _PM.variable_branch_power(pm)

    objective_min_fuel_and_consensus!(pm, objective_admm_coordinated!)

end


"ADMM algorithm objective"
function objective_admm_coordinated!(pm::AbstractPowerModel)

    ## ADMM parameters
    alpha = pm.data["alpha"]

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
                v_central = primal_variable[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]
 
                objective += alpha/2 * (v - v_central)^2 + v_dual * (v - v_central)
            end
        end
    end

    JuMP.@objective(pm.model, Min,  objective)
end



function update_admm_coordinated!(data::Dict{String, <:Any})

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
                v_local =  primal_variable[area][variable][idx]
                v_dual = dual_variable[area][variable][idx]

                data["shared_dual"][area][variable][idx]= v_dual  + alpha * (v_primal - v_local)
            end
        end
    end

end