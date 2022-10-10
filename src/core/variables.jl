###############################################################################
#   Variable initilization and updating for all distirbuted OPF algorithms    #
###############################################################################

# Fully distributed initialization
"initialize the primal shared variables"
function initialize_shared_variable!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id

    data["shared_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", method)

    data["received_shared_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", method)
end

"initialize the dual shared variables"
function initialize_dual_variable!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id =  get_areas_id(data)
   
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id
    data["dual_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "dual_variable", method)

end


# Coordinated distributed initialization
"initialize the primal shared variables for coordinator"
function initialize_shared_variable_coordinator!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id = get_areas_id(data)

    data["shared_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", method)

    data["received_shared_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", method)


    # data["shared_variable"] = _initialize_shared_variable(data, model_type, area_id ,[areas_id; area_id], "shared_variable", method)

end

"initialize the dual shared variables for coordinator"
function initialize_dual_variable_coordinator!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
   
    data["dual_variable"] = _initialize_shared_variable(data, model_type, area_id, areas_id, "dual_variable", method)
    
end

# Local systems initialization
"initialize the primal shared variables"
function initialize_shared_variable_local!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)
    
    data["shared_variable"] = _initialize_shared_variable(data, model_type, area_id, [0], "shared_variable", method)

    data["received_shared_variable"] = _initialize_shared_variable(data, model_type, area_id, [0], "shared_variable", method)

    # data["shared_variable"] = _initialize_shared_variable(data, model_type, area_id ,[0, area_id], "shared_variable", method)

end

"initialize the dual shared variables"
function initialize_dual_variable_local!(data::Dict{String, <:Any}, model_type::DataType; method::String="flat")
    area_id = get_area_id(data)

    data["dual_variable"] = _initialize_shared_variable(data, model_type, area_id ,[0], "dual_variable", method)
end

# Template for variable shared
"initialize shared variable dictionary"
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
function initial_value(data::Dict{String, <:Any}, variable::String, idx::Int64, area::Int64, method::String="flat", shared_variable::String="shared_variable")::Float64

    if method in ["flat" , "flat_start"]
        if shared_variable == "shared_variable" && var in ["vm", "w", "wr"]
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

"""
    initialize_all_variable(data::Dict{String, <:Any}, model_type::DataType)

return a dictionary contains all the problem variables. Can be used to store the last solution

# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
"""
function initialize_all_variable(data::Dict{String, <:Any}, model_type::DataType)
    bus_variables_name, branch_variables_name, gen_variables_name = variable_names(model_type)
    all_variables = Dict{String, Dict}()
    for variable in bus_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, parse(Int64,idx), 0, "flat", "shared_variable")  for idx in keys(data["bus"])])
    end
    for variable in branch_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, parse(Int64,idx), 0, "flat", "shared_variable")  for idx in keys(data["branch"])])
    end
    for variable in gen_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, parse(Int64,idx), 0, "flat", "shared_variable")  for idx in keys(data["gen"])])
    end
    return all_variables
end

function initialize_solution!(data::Dict{String, <:Any}, model_type::DataType)
    data["solution"] = initialize_all_variable(data, model_type)
end

"return JuMP variable object from PowerModel object"
function _var(pm::AbstractPowerModel, key::String, idx::String)
    bus_variables_name, branch_variables_name, gen_variables_name = variable_names(typeof(pm))
    idx = parse(Int64,idx)
    if key in bus_variables_name || key in gen_variables_name
        var = _PM.var(pm, Symbol(key), idx)
    elseif key in branch_variables_name
        branch = _PM.ref(pm, :branch, idx)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]

        if key in ["pf", "qf"]
            var = _PM.var(pm, Symbol(key[1]),  (idx, f_bus, t_bus))
        elseif key in ["pt", "qt"]
            var = _PM.var(pm, Symbol(key[1]),  (idx, t_bus, f_bus))
        else
            var = _PM.var(pm, Symbol(key), (f_bus, t_bus))
        end
    end

    return var
end

"idinifiy the shared bus and branch variables names"
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

"idinifiy all the variables names"
function variable_names(model_type::DataType)
    if model_type <: Union{DCPPowerModel, DCMPPowerModel}
        return ["va"], ["pf"], ["pg"]
    elseif model_type <: NFAPowerModel
        return [], ["pf"], ["pg"]
    elseif model_type <: DCPLLPowerModel
        return ["va"], ["pf", "pt"], ["pg"]
    elseif model_type <: ACPPowerModel
        return ["va", "vm"], ["pf", "pt", "qf", "qt"], ["pg", "qg"]
    elseif model_type <: ACRPowerModel
        return ["vr", "vi"], ["pf", "pt", "qf", "qt"], ["pg", "qg"]
    elseif model_type <: ACTPowerModel
        return ["w", "va"], ["pf", "pt", "qf", "qt", "wr", "wi"], ["pg", "qg"]
    elseif model_type <: Union{SOCWRPowerModel, SOCWRConicPowerModel, SDPWRMPowerModel, SparseSDPWRMPowerModel }
        return ["w"], ["pf", "pt", "qf", "qt", "wr", "wi"], ["pg", "qg"]
    elseif model_type <: QCRMPowerModel
        return ["vm", "va" , "w"], ["pf", "pt", "qf", "qt", "wr", "wi", "vv", "ccm", "cs", "si", "td"], ["pg", "qg"]
    elseif model_type <: AbstractPowerModel
        error("PowerModel type is not supported yet!")
    else
        error("model_type $model_type is not PowerModel type!")
    end
end