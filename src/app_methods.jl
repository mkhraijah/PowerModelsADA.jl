###############################################################################
#                     Build methods for APP algorithm                        #
###############################################################################



## build method for Distributed PowerModel using ADMM algorithm
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

    update_dual_app!(pm)

    objective_min_fuel_and_consensus!(pm, objective_app!)
end

## method to set the APP algorithm objective
function objective_app!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.setting["alpha"]
    # use beta if defined in setting or use 2α
    if haskey(pm.setting, "beta")
        beta = pm.setting["beta"]
    else
        beta = 2*alpha
    end
    # use gamma if defined in setting or use α
    if haskey(pm.setting, "gamma")
        gamma = pm.setting["gamma"]
    else
        gamma = alpha
    end

    ## data
    area = get_area_id(pm)
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## objective function
    objective = JuMP.objective_function(pm.model) + sum(beta/2 * (var(pm, j, k) - primal_variable[area][j][k])^2 + var(pm, j, k) * dual_variable[i][j][k] + gamma * var(pm, j, k) * (primal_variable[area][j][k] - primal_variable[i][j][k]) for i in keys(primal_variable) if i != area for j in keys(primal_variable[i]) for k in keys(primal_variable[i][j]))

    JuMP.@objective(pm.model, Min,  objective)
end


## method to update the dual variable value

function update_dual_app!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.setting["alpha"]

    ## data
    area_id = get_area_id(pm)
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## update dual variable
    for i in keys(dual_variable)
        for j in keys(dual_variable[i])
            for k in keys(dual_variable[i][j])
                dual_variable[i][j][k] = dual_variable[i][j][k] + alpha * (primal_variable[area_id][j][k] - primal_variable[i][j][k])
            end
        end
    end

end

function run_dopf_app(data, pf_model, optimizer, alpha::Real=1000; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true)
    run_dopf(data, pf_model, build_dopf_app, optimizer, alpha, tol = tol, max_iteration=max_iteration,verbose =verbose)
end
