## Method to handle shared variables storing and updating
## shared variables (primal and dual) are stored in pm.ext
function variable_shared(pm::AbstractPowerModel)
    primal_variable_shared!(pm)
    dual_variable_shared!(pm)
end

function primal_variable_shared!(pm::AbstractPowerModel)

    area_id = get_area_id(pm)
    areas_id = get_areas_id(pm)
    nodal_variables_name, cross_variables_name = variable_shared_names(pm)
    shared_bus, shared_branch = get_shared_component(pm, area_id)

    pm.ext[:primal_shared_variable] = Dict{Int64, Any}()
    ## Nodal varaibles
    for i in areas_id
        pm.ext[:primal_shared_variable][i] =  Dict{Symbol, Any}()
        for k in nodal_variables_name
            pm.ext[:primal_shared_variable][i][k] = Dict{Int64, Float64}([j=> 0 for j in shared_bus[i]])
        end

    ## Cross variables
        for k in cross_variables_name
            if k in [:p, :q]
                pm.ext[:primal_shared_variable][i][k] = merge( Dict{Tuple{Int64, Int64, Int64}, Float64}([(j,pm.data["branch"]["$j"]["f_bus"],pm.data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if pm.data["branch"]["$j"]["br_status"] ==1 ]), Dict{Tuple{Int64, Int64, Int64}, Float64}([(j,pm.data["branch"]["$j"]["t_bus"],pm.data["branch"]["$j"]["f_bus"]) => 0 for j in shared_branch[i] if pm.data["branch"]["$j"]["br_status"] ==1 ]))
            else
                pm.ext[:primal_shared_variable][i][k] = Dict{Tuple{Int64, Int64}, Any}([ (pm.data["branch"]["$j"]["f_bus"],pm.data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if (pm.data["branch"]["$j"]["br_status"] == 1)])
            end
        end
    end
end

## primal shared variable for DC case
function primal_variable_shared!(pm::AbstractDCPModel)
    area_id = get_area_id(pm)
    areas_id = get_areas_id(pm)
    nodal_variables_name, cross_variables_name = variable_shared_names(pm)
    shared_bus, shared_branch = get_shared_component(pm, area_id)

    pm.ext[:primal_shared_variable] = Dict{Int64, Any}()
    ## Nodal varaibles
    for i in areas_id
        pm.ext[:primal_shared_variable][i] =  Dict{Symbol, Any}()
        for k in nodal_variables_name
            pm.ext[:primal_shared_variable][i][k] = Dict{Int64, Float64}([j=> 0 for j in shared_bus[i]])
        end

    ## Cross variables
        for k in cross_variables_name
            pm.ext[:primal_shared_variable][i][k] = Dict{Tuple{Int64, Int64, Int64}, Float64}([(j,pm.data["branch"]["$j"]["f_bus"],pm.data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if pm.data["branch"]["$j"]["br_status"] == 1 ])
        end
    end

end

function dual_variable_shared!(pm::AbstractPowerModel)
    area_id = get_area_id(pm)
    areas_id = get_areas_id(pm)

    pm.ext[:dual_shared_variable] = Dict{Int64, Any}()
    for i in areas_id
        if i != area_id
            pm.ext[:dual_shared_variable][i] = deepcopy(pm.ext[:primal_shared_variable][i])
        end
    end
end

## method to update pm.ext[:primal_shared_variable] after obtaining a solution at each iteraton
function update_primal_variable!(pm::AbstractPowerModel)
    area_id = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    for i in keys(primal_variable[area_id])
        for j in keys(primal_variable[area_id][i])
            primal_variable[area_id][i][j] = JuMP.value(var(pm,i,j))
        end
    end
end

""" is there a better way to handle multiple algorithms instead of if statement? """
function update_dual_variable!(pm::AbstractPowerModel)

    Algorithm = pm.setting["distributed_algorithm"]
    if Algorithm == "APP"
        update_dual_APP!(pm)
    elseif Algorithm == "ADMM"
        update_dual_ADMM!(pm)
    elseif Algorithm == "ATC"
        update_dual_ATC!(pm)
    end
end

function update_dual_APP!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.ext[:parameter][1]

    ## Data
    area_id = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## Update Lagrange Multipliers
    for i in keys(dual_variable)
        for j in keys(dual_variable[i])
            for k in keys(dual_variable[i][j])
                dual_variable[i][j][k] = dual_variable[i][j][k] + alpha * (primal_variable[area_id][j][k] - primal_variable[i][j][k])
            end
        end
    end
end


function update_dual_ADMM!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.ext[:parameter][1]

    ## Data
    area_id = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## Update Lagrange Multipliers
    for i in keys(dual_variable)
        for j in keys(dual_variable[i])
            for k in keys(dual_variable[i][j])
                dual_variable[i][j][k] = dual_variable[i][j][k] + alpha * (primal_variable[area_id][j][k] - (primal_variable[area_id][j][k] +primal_variable[i][j][k])/2 )
            end
        end
    end
end


function update_dual_ATC!(pm::AbstractPowerModel)

    ## APP parameters
    beta = pm.ext[:beta]

    ## Data
    area_id = get_area_id(pm)
    primal_variable = pm.ext[:primal_shared_variable]
    dual_variable = pm.ext[:dual_shared_variable]

    ## Update Lagrange Multipliers
    for i in keys(dual_variable)
        for j in keys(dual_variable[i])
            for k in keys(dual_variable[i][j])
                dual_variable[i][j][k] = dual_variable[i][j][k] + 2 * beta^2 * (primal_variable[area_id][j][k] - (primal_variable[area_id][j][k]+primal_variable[i][j][k])/2 )
            end
        end
    end

    ## Update ATC parameter
    pm.ext[:beta] *= pm.ext[:parameter][1]

end

## Methods to idinifiy the nodal and cross variables names
variable_shared_names(pm::AbstractDCPModel) = [:va], [:p]

variable_shared_names(pm::AbstractACPModel) = [:va, :vm], [:p, :q]

variable_shared_names(pm::AbstractACRModel) = [:vr, :vi], [:p, :q]

variable_shared_names(pm::AbstractSOCWRModel) = [:w], [:p, :q, :wr, :wi]

variable_shared_names(pm::AbstractQCRMPowerModel) = [:vm, :va , :w], [:p, :q, :wr, :wi, :vv, :ccm, :cs, :si, :td]

variable_shared_names(pm::AbstractSDPWRMModel) = [:w], [:p, :q, :wr, :wi]
