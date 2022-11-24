###############################################################################
#             OPF problem variable, objetive, and constraints                 #
###############################################################################

"define OPF problem variable"
function variable_opf(pm::AbstractPowerModel)
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)
end

"define objective function using PowerModels and algorithm-specific objective"
function objective_min_fuel_and_consensus!(pm::AbstractPowerModel, objective_method::Function=no_objective)
    # if subsystem has generator minimize the cost of generator and consistency otherwise minimize consistency only
    if isempty(pm.data["gen"])
        objective = objective_method(pm)
    else
        _PM.objective_min_fuel_and_flow_cost(pm) 
        objective = JuMP.objective_function(pm.model) + objective_method(pm)
    end
    
    JuMP.@objective(pm.model, Min, objective)
end

"no objective function case"
function no_objective(pm::AbstractPowerModel)
    # do nothing
end


"define OPF problem constraints"
function constraint_opf(pm::AbstractPowerModel)
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
end
