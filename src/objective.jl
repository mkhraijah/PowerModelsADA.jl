## methods to handle objective function update at every iteration

## method to update the objective function at each iteration
""" will be changed to instantiate the subsystem at each iteration using the buidling function """

function update_objective!(pm::AbstractPowerModel; instantiate::Bool=false)

    if isa(pm, _PM.AbstractConicModel) ||  _PM.check_cost_models(pm) == 1
        instantiate = true
    end

    if instantiate
        pf_model = typeof(pm)
        pm_tmp = _PM.instantiate_model(pm.data, pf_model, build_dopf)
        update_pm!(pm,pm_tmp) ## retrive the previosully stored information
        optimizer = pm.ext[:optimizer]
        JuMP.set_optimizer(pm.model, optimizer)
        JuMP.set_silent(pm.model)
    end

    objective_min_fuel_and_consensus!(pm)
end

function update_pm!(pm, pm_temp)
    pm.model = pm_temp.model
    pm.var = pm_temp.var
    pm.sol = pm_temp.sol
end

## wrapper for objective function definition from PowerModels
function objective_min_fuel_and_consensus!(pm::AbstractPowerModel)

    # if subsystem has generator minimize the cost of generator
    if !isempty(pm.data["gen"])
        _PM.objective_min_fuel_and_flow_cost(pm)
    else
        JuMP.set_objective_function(pm.model, 0)
    end
    # set consensus penality based on distributed algorithm
    objective_consensus!(pm)

end

## consensus objective function update based on the selected algorithm
function objective_consensus!(pm::AbstractPowerModel)
    Algorithm = pm.setting["distributed_algorithm"]
    if Algorithm == "APP"
        objective_consensus_APP!(pm)
    elseif Algorithm == "ADMM"
        objective_consensus_ADMM!(pm)
    elseif Algorithm == "ATC"
        objective_consensus_ATC!(pm)
    end
end

##
function objective_consensus_APP!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.ext[:parameter][1]
    if length(pm.ext[:parameter]) > 2
        beta = pm.ext[:parameter][2] # use beta if defined in setting or use 2α
        gamma = pm.ext[:parameter][3] # use gamma if defined in setting or use α
    else
        beta = 2*alpha
        gamma = alpha
    end
    ## data
    area = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## objective function
    objective = JuMP.objective_function(pm.model) + sum(beta/2 * (var(pm, j, k) - primal_variable[area][j][k])^2 + var(pm, j, k) * dual_variable[i][j][k] + gamma * var(pm, j, k) * (primal_variable[area][j][k] - primal_variable[i][j][k]) for i in keys(primal_variable) if i != area for j in keys(primal_variable[i]) for k in keys(primal_variable[i][j]))

    JuMP.@objective(pm.model, Min,  objective)
end

function objective_consensus_ADMM!(pm::AbstractPowerModel)
    ## ADMM parameters
    alpha = pm.ext[:parameter][1]

    ## Data
    area = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## objective function
    objective = JuMP.objective_function(pm.model) + sum(dual_variable[i][j][k] * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2) + alpha/2 * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2)^2 for i in keys(primal_variable) if i != area for j in keys(primal_variable[i]) for k in keys(primal_variable[i][j]))

    JuMP.@objective(pm.model, Min,  objective)
end


function objective_consensus_ATC!(pm::AbstractPowerModel)
    ## ATC parameters
    beta = pm.ext[:beta]
    ## Data
    area = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## objective function
    objective = JuMP.objective_function(pm.model) + sum(dual_variable[i][j][k] * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2) + (beta * (var(pm, j, k) - (primal_variable[area][j][k] + primal_variable[i][j][k])/2))^2 for i in keys(primal_variable) if i != area for j in keys(primal_variable[i]) for k in keys(primal_variable[i][j]))

    JuMP.@objective(pm.model, Min,  objective)
end
