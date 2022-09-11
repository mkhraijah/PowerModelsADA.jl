###############################################################################
#   Variable initilization and updating for all distirbuted OPF algorithms    #
###############################################################################

"initialize the shared variables"
function initialize_variable_shared!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    initialize_primal_variable_shared!(data, model_type)
    initialize_dual_variable_shared!(data, model_type)
end

"initialize the primal shared variables"
function initialize_primal_variable_shared!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
   
    data["shared_primal"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_primal", method)
end

"initialize the dual shared variables"
function initialize_dual_variable_shared!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
   
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id
    data["shared_dual"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_dual", method)

end

"initialize the shared variables for coordinator"
function initialize_variable_shared_coordinator!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    initialize_primal_variable_shared_coordinator!(data, model_type)
    initialize_dual_variable_shared_coordinator!(data, model_type)
end

"initialize the primal shared variables for coordinator"
function initialize_primal_variable_shared_coordinator!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)

    data["shared_primal"] = _initialize_shared_variable(data, model_type, area_id ,[areas_id; area_id], "shared_primal", method)

end

"initialize the dual shared variables for coordinator"
function initialize_dual_variable_shared_coordinator!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
   
    data["shared_dual"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_dual", method)
    
end

"initialize the shared variables"
function initialize_variable_shared_local!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    initialize_primal_variable_shared_local!(data, model_type)
    initialize_dual_variable_shared_local!(data, model_type)
end

"initialize the primal shared variables"
function initialize_primal_variable_shared_local!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    
    data["shared_primal"] = _initialize_shared_variable(data, model_type, area_id ,[0, area_id], "shared_primal", method)

end

"initialize the dual shared variables"
function initialize_dual_variable_shared_local!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)

    data["shared_dual"] = _initialize_shared_variable(data, model_type, area_id ,[0], "shared_dual", method)
end

"initlize shared variable dictionary"
function _initialize_shared_variable(data::Dict{String, <:Any}, model_type::DataType, from::Int64 ,to::Vector{Int64}, shared_variable::String, method::String="flat")
    bus_variables_name, branch_variables_name =  variable_shared_names(model_type)
    shared_bus, shared_branch =  get_shared_component(data, from)

    return Dict{String, Any}([
        string(area) => Dict{String, Any}(
            vcat(
                [variable => Dict{String, Any}([string(idx) => initial_value(data, variable, idx, area, method, shared_variable) for idx in shared_bus[area]]) for variable in bus_variables_name],
                [variable => Dict{String, Any}([string(idx) => initial_value(data, variable, idx, area, method, shared_variable) for idx in shared_branch[area]]) for variable in branch_variables_name]
            )
        )
    for area in to])
end

"""
    initial_value(data::Dict{String, <:Any}, var::String, idx::Int, method::String="flat")

assign initial value based on initialization method

# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- var::String : variable names
- idx::Int : variable index
- method::String : initialization method ("flat", "previous_solution")
"""
function initial_value(data::Dict{String, <:Any}, variable::String, idx::Int64, area::Int64, method::String="flat", shared_variable::String="shared_primal")

    if method in ["flat" , "flat_start"]
        if shared_variable == "shared_primal" && var in ["vm", "w", "wr"]
            return 1
        else
            return 0
        end
    elseif method in ["previous", "previous_solution", "warm", "warm_start"]
        if !haskey(data, shared_variable)
            error("no previous solution exist to use warm start")
        end
        return data[shared_variable][string(area)][variable][idx]
    else
        error("the initlization method is not supported")
    end
end


"idinifiy the nodal and cross variables names"
function variable_shared_names(model_type::DataType)
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
    else
        error("PowerModel type is not supported yet!")
    end
end