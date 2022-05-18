###############################################################################
#                     Build methods for ATC algorithm                        #
###############################################################################



## build method for Distributed PowerModel using ADMM algorithm
function build_dopf_atc(pm::AbstractPowerModel)

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

    update_dual_atc!(pm)

    objective_min_fuel_and_consensus!(pm, objective_atc!)
end

## method to set the ATC algorithm objective
function objective_atc!(pm::AbstractPowerModel)

    ## ATC parameters
    beta = pm.data["beta"]

    ## data
    area = get_area_id(pm)
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## objective function
    objective = JuMP.objective_function(pm.model) + sum(dual_variable[i][j][k] * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2) + (beta * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2))^2 for i in keys(primal_variable) if i != area for j in keys(primal_variable[i]) for k in keys(primal_variable[i][j]))

    JuMP.@objective(pm.model, Min,  objective)
end


## method to update the dual variable value
function update_dual_atc!(pm::AbstractPowerModel)

    ## ATC parameters
    if pm.data["iteration"] == 1
        pm.data["beta"] = 1
    end

    alpha = pm.setting["alpha"]
    beta = pm.data["beta"]

    ## data
    area_id = get_area_id(pm)
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## update dual variable
    for i in keys(dual_variable)
        for j in keys(dual_variable[i])
            for k in keys(dual_variable[i][j])
                dual_variable[i][j][k] = dual_variable[i][j][k] + 2 * beta^2 * (primal_variable[area_id][j][k] - (primal_variable[area_id][j][k]+primal_variable[i][j][k])/2 )
            end
        end
    end

    ## update ATC parameter
    pm.data["beta"] *= alpha

end


function run_dopf_atc(data, pf_model, optimizer, alpha::Real=1.05; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true)
    run_dopf(data, pf_model, build_dopf_atc, optimizer, alpha, tol = tol, max_iteration=max_iteration,verbose =verbose)
end
