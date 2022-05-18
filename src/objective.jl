###############################################################################
#           Objective update for all distirbuted OPF algorithms               #
###############################################################################

## wrapper for objective function definition from PowerModels and algorithm specific objective
function objective_min_fuel_and_consensus!(pm::AbstractPowerModel, objective_method::Function)

    # if subsystem has generator minimize the cost of generator
    if !isempty(pm.data["gen"])
        _PM.objective_min_fuel_and_flow_cost(pm)
    else
        JuMP.set_objective_function(pm.model, 0)
    end

    # set consensus penality based on distributed algorithm
    objective_method(pm)

end
