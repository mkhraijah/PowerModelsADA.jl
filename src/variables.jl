###############################################################################
#   Variable initilization and updating for all distirbuted OPF algorithms    #
###############################################################################

## initialize the shared variables
function variable_shared(data::Dict{String, <:Any}, pf_model)
    primal_variable_shared!(data, pf_model)
    dual_variable_shared!(data)
end

## initialize the primal shared variables
function primal_variable_shared!(data::Dict{String, <:Any}, pf_model)

    ## Check the power flow model
    pf_model= pf_formulation(pf_model)

    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
    nodal_variables_name, cross_variables_name = variable_shared_names(pf_model)
    shared_bus, shared_branch = get_shared_component(data, area_id)

    data["shared_primal"] = Dict{Int64, Any}()

    for i in areas_id
        data["shared_primal"][i] =  Dict{Symbol, Any}()
        ## Nodal varaibles
        for k in nodal_variables_name
            data["shared_primal"][i][k] = Dict{Int64, Float64}([j=> 0 for j in shared_bus[i] if data["bus"]["$j"]["bus_type"] != 4])
        end

        ## Cross variables
        if pf_model <: AbstractDCPModel
            for k in cross_variables_name
                data["shared_primal"][i][k] = Dict{Tuple{Int64, Int64, Int64}, Float64}([(j, data["branch"]["$j"]["f_bus"], data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if data["branch"]["$j"]["br_status"] == 1 ])
            end
        else
            for k in cross_variables_name
                if k in [:p, :q]
                    data["shared_primal"][i][k] = merge( Dict{Tuple{Int64, Int64, Int64}, Float64}([(j, data["branch"]["$j"]["f_bus"], data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if data["branch"]["$j"]["br_status"] == 1 ]), Dict{Tuple{Int64, Int64, Int64}, Float64}([(j, data["branch"]["$j"]["t_bus"], data["branch"]["$j"]["f_bus"]) => 0 for j in shared_branch[i] if data["branch"]["$j"]["br_status"] == 1 ]))
                else
                    data["shared_primal"][i][k] = Dict{Tuple{Int64, Int64}, Any}([ (data["branch"]["$j"]["f_bus"], data["branch"]["$j"]["t_bus"]) => 0 for j in shared_branch[i] if (data["branch"]["$j"]["br_status"] == 1)])
                end
            end
        end
    end
end

## initialize the dual shared variables
function dual_variable_shared!(data::Dict{String, <:Any})
    area_id = get_area_id(data)
    areas_id = get_areas_id(data)

    data["shared_dual"] = Dict{Int64, Any}()
    for i in areas_id
        if i != area_id
            data["shared_dual"][i] = deepcopy(data["shared_primal"][i])
        end
    end
end


## method to update primal variables after obtaining a solution at each iteraton
function update_primal_variable!(pm::AbstractPowerModel)
    area_id = get_area_id(pm)
    primal_variable = pm.data["shared_primal"]
    for i in keys(primal_variable[area_id])
        for j in keys(primal_variable[area_id][i])
            primal_variable[area_id][i][j] = JuMP.value(var(pm,i,j))
        end
    end
end


## Methods to idinifiy the nodal and cross variables names
function variable_shared_names(pf_model)
    if pf_model <: AbstractDCPModel
        return [:va], [:p]
    elseif pf_model <: AbstractACPModel
        return [:va, :vm], [:p, :q]
    elseif pf_model <: AbstractACRModel
        return [:vr, :vi], [:p, :q]
    elseif pf_model <: AbstractSOCWRModel
        return [:w], [:p, :q, :wr, :wi]
    elseif pf_model <: AbstractQCRMPowerModel
        return [:vm, :va , :w], [:p, :q, :wr, :wi, :vv, :ccm, :cs, :si, :td]
    elseif pf_model <: AbstractSDPWRMModel
        return [:w], [:p, :q, :wr, :wi]
    end
end
