###############################################################################
#   Variable initilization and updating for all distirbuted OPF algorithms    #
###############################################################################

"initialize the shared variables"
function variable_shared(data::Dict{String, <:Any}, model_type)
    primal_variable_shared!(data, model_type)
    dual_variable_shared!(data, model_type)
end

"initialize the primal shared variables"
function primal_variable_shared!(data::Dict{String, <:Any}, model_type)
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
    bus_variables_name, branch_variables_name =  variable_shared_names(model_type)
    shared_bus, shared_branch =  get_shared_component(data, area_id)

    # data["shared_primal"] = Dict{String, Any}()

    data["shared_primal"] = Dict{String, Any}(
        Dict{String, Any}([
            string(area) => Dict{String, Any}(
                vcat(
                    [variable => Dict{String, Any}([string(idx) => 0 for idx in shared_bus[area]]) for variable in bus_variables_name],
                    [variable => Dict{String, Any}([string(idx) => 0 for idx in shared_branch[area]]) for variable in branch_variables_name]
                )
            )
        for area in areas_id])
    )
end

"initialize the dual shared variables"
function dual_variable_shared!(data::Dict{String, <:Any}, model_type)
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
    bus_variables_name, branch_variables_name =  variable_shared_names(model_type)
    shared_bus, shared_branch =  get_shared_component(data, area_id)

    data["shared_dual"] = Dict{String, Any}(
        Dict{String, Any}([
            string(area) => Dict{String, Any}(
                vcat(
                    [variable => Dict{String, Any}([string(idx) => 0 for idx in shared_bus[area]]) for variable in bus_variables_name],
                    [variable => Dict{String, Any}([string(idx) => 0 for idx in shared_branch[area]]) for variable in branch_variables_name]
                )
            )
        for area in areas_id if area != area_id])
    )

end

"idinifiy the nodal and cross variables names"
function variable_shared_names(model_type)
    if model_type <: Union{DCPPowerModel, DCMPPowerModel}
        return ["va"], ["pf"]
    elseif model_type <: NFAPowerModel
        return [], ["pf"]
    elseif model_type <: DCPLLPowerModel
        return ["va"], ["pf", "pt"]
    elseif model_type <: ACPPowerModel
        return ["va", "vm"], ["pf", "pt", "qf", "qt"]
    elseif model_type <: ACRPowerModel
        return ["vr", "vi"], ["pf", "pt", "qf", "qt"]
    elseif model_type <: ACTPowerModel
        return ["w", "va"], ["pf", "pt", "qf", "qt", "wr", "wi"]
    elseif model_type <: Union{SOCWRPowerModel, SOCWRConicPowerModel, SDPWRMPowerModel, SparseSDPWRMPowerModel }
        return ["w"], ["pf", "pt", "qf", "qt", "wr", "wi"]
    elseif model_type <: QCRMPowerModel
        return ["vm", "va" , "w"], ["pf", "pt", "qf", "qt", "wr", "wi", "vv", "ccm", "cs", "si", "td"]
    end
end
