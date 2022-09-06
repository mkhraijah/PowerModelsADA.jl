###############################################################################
#           Objective update for all distirbuted OPF algorithms               #
###############################################################################

"define objective function from PowerModels and algorithm specific objective"
function objective_min_fuel_and_consensus!(pm::AbstractPowerModel, objective_method::Function=no_objective)

    # if subsystem has generator minimize the cost of generator
    if !isempty(pm.data["gen"])
        _PM.objective_min_fuel_and_flow_cost(pm)
    end
    objective = JuMP.objective_function(pm.model)

    # add consensus penality based on distributed algorithm
    objective += objective_method(pm)
    JuMP.@objective(pm.model, Min,  objective)

end

function no_objective(pm::AbstractPowerModel)
    # do nothing
end