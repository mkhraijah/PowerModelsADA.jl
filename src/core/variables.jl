###############################################################################
#   Variable initilization and updating for all distirbuted OPF algorithms    #
###############################################################################

## initialize the shared variables
function variable_shared(data::Dict{String, <:Any}, model_type)
    primal_variable_shared!(data, model_type)
    dual_variable_shared!(data)
end

## initialize the primal shared variables
function primal_variable_shared!(data::Dict{String, <:Any}, model_type)

    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
    bus_variables_name, branch_variables_name =  variable_shared_names(model_type)
    shared_bus, shared_branch =  get_shared_component(data, area_id)

    data["shared_primal"] = Dict{String, Any}()

    for i in areas_id
        data["shared_primal"][string(i)] =  Dict{String, Any}(["bus" => Dict{String, Any}(), "branch" => Dict{String, Any}()])

        ## bus varaibles
        for j in shared_bus[i]
            if data["bus"][string(j)]["bus_type"] != 4
                data["shared_primal"][string(i)]["bus"][string(j)] = Dict{String, Any}()
                for k in bus_variables_name
                    data["shared_primal"][string(i)]["bus"][string(j)][k] = 0
                end
            end
        end

        ## branch variables
        for j in shared_branch[i]
            if data["branch"][string(j)]["br_status"] == 1
                data["shared_primal"][string(i)]["branch"][string(j)] = Dict{String, Any}()
                for k in branch_variables_name
                    data["shared_primal"][string(i)]["branch"][string(j)][k] = 0
                end
            end
        end
    end
end

## initialize the dual shared variables
function dual_variable_shared!(data::Dict{String, <:Any})
    area_id = get_area_id(data)
    areas_id = get_areas_id(data)

    data["shared_dual"] = Dict{String, Any}()
    for i in areas_id
        if i != area_id
            data["shared_dual"][string(i)] = deepcopy(data["shared_primal"][string(i)])
        end
    end
end


## method to update primal variables after obtaining a solution at each iteraton
function update_shared_primal!(data::Dict{String,<:Any}, solution::Dict{String,<:Any})
    area_id = string(get_area_id(data))
    shared_primal = data["shared_primal"][area_id]
    _update_shared_primal!(shared_primal, solution)
end

function _update_shared_primal!(data::Dict{String,<:Any}, new_data::Dict{String,<:Any})

    for (key, new_v) in new_data
        if haskey(data, key)
            v = data[key]
            if isa(v, Dict) && isa(new_v, Dict)
                _update_shared_primal!(v, new_v)
            else
                data[key] = new_v
            end
        end
    end
end

## Methods to idinifiy the nodal and cross variables names
function variable_shared_names(model_type)
    if model_type <: DCPPowerModel
        return ["va"], ["pf"]
    elseif model_type <: ACPPowerModel
        return ["va", "vm"], ["pf", "pt", "qf", "qt"]
    elseif model_type <: ACRPowerModel
        return ["vr", "vi"], ["pf", "pt", "qf", "qt"]
    elseif model_type <: SOCWRPowerModel
        return ["w"], ["pf", "pt", "qf", "qt", "wr", "wi"]
    elseif model_type <: QCRMPowerModel
        return ["vm", "va" , "w"], ["pf", "pt", "qf", "qt", "wr", "wi", "vv", "ccm", "cs", "si", "td"]
    elseif model_type <: SDPWRMPowerModel
        return ["w"], ["pf", "pt", "qf", "qt", "wr", "wi"]
    end
end
