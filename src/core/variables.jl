###############################################################################
#   Variable initialization and updating for all distributed OPF algorithms   #
###############################################################################

# Template for variable shared
"initialize shared variable dictionary"
function initialize_shared_variable(data::Dict{String, <:Any}, model_type::DataType, from::Int64 ,to::Vector{Int64}, dics_name::String="shared_variable", initialization_method::String="flat", value::Float64=0.0)
    bus_variables_name, branch_variables_name = variable_shared_names(model_type)
    shared_bus, shared_branch = get_shared_component(data, from)

    if initialization_method in ["previous", "previous_solution", "warm", "warm_start"]
        if !haskey(data, dics_name)
            error("no previous solutions exist to use warm start")
        else
            variables_dics = data[dics_name]
        end
    else
        variables_dics = Dict{String, Any}([
            string(area) => Dict{String, Any}(
                vcat(
                    [variable => Dict{String, Any}([string(idx) => initial_value(data, variable, string(idx), initialization_method, value) for idx in shared_bus[area]]) for variable in bus_variables_name],
                    [variable => Dict{String, Any}([string(idx) => initial_value(data, variable, string(idx), initialization_method, value) for idx in shared_branch[area]]) for variable in branch_variables_name]
                )
            )
        for area in to])
    end

    return variables_dics
end

function initialize_shared_variable(data::Dict{String, <:Any}, model_type::DataType, from::Int64 ,to::Int64, dics_name::String="shared_variable", initialization_method::String="flat")
    initialize_shared_variable(data, model_type, from, [to], dics_name, initialization_method)
end

"""
    initial_value(variable::String, initialization_method::String="flat")

assign initial value based on initialization method

# Arguments:
- variable::String : variable names
- initialization_method::String="flat : ("flat", "previous_solution")
"""
function initial_value(data::Dict{String, <:Any}, variable::String, idx::String, initialization_method::String="flat", value::Float64=0.0)::Float64
    if initialization_method in ["previous", "previous_solution", "warm", "warm_start"]
        return previous_value(data, variable, idx)
    elseif initialization_method in ["flat" , "flat_start"]
        return initial_value(variable)
    elseif initialization_method in ["constant"]
        return value
    else
        return 0.0
    end
end

function initial_value(variable::String)::Float64
    if variable in ["vm", "w", "wr"]
        return 1.0
    else
        return 0.0
    end
end

function previous_value(data::Dict{String, <:Any}, variable::String, idx::String)::Float64
    
    if variable in ["vm", "va"]
        return data["bus"][idx][variable]
    elseif variable in ["w"]
        if haskey(data["bus"][idx], "w")
            return data["bus"][idx]["w"]
        else
            return data["bus"][idx]["vm"]^2
        end
    elseif variable in ["pf", "pt", "qf", "qt","wr", "wi", "vv", "ccm", "cs", "si", "td"]
        if haskey(data["bus"][idx], variable)
            return data["branch"][idx][variable]
        else
            error("no previous solutions exist to use warm start or the PowerModel is not supported")
        end
    elseif ["pg", "qg"]
        return data["gen"][idx][variable]
    else
        error("no previous solutions exist to use warm start or the PowerModel is not supported")
    end

end

"""
    initialize_all_variable(data::Dict{String, <:Any}, model_type::DataType, dics_name::String="solution", initialization_method::String="flat")

return a dictionary contains all the problem variables. can be used to store the solutions.

# Arguments:
- data::Dict{String, <:Any} : area data
- model_type::DataType : power flow formulation (PowerModel type)
- dics_name::String="solution" : location of existing dicrionary to be used to worm start the output
- initialization_method::String="flat" : "flat" or "worm" initialization
"""
function initialize_all_variable(data::Dict{String, <:Any}, model_type::DataType, initialization_method::String="flat")
    bus_variables_name, branch_variables_name, gen_variables_name = variable_names(model_type)

    all_variables = Dict{String, Dict}()
    for variable in bus_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, idx, initialization_method) for idx in keys(data["bus"])])
    end
    for variable in branch_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, idx, initialization_method) for idx in keys(data["branch"])])
    end
    for variable in gen_variables_name
        all_variables[variable] = Dict([idx => initial_value(data, variable, idx, initialization_method) for idx in keys(data["gen"])])
    end

    return all_variables
end

function initialize_solution!(data::Dict{String, <:Any}, model_type::DataType, initialization_method::String="flat")
    data["solution"] = initialize_all_variable(data, model_type, initialization_method)
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

"identifythe shared bus and branch variables names"
function variable_shared_names(model_type::DataType)
    if model_type <: Union{DCPPowerModel, DCMPPowerModel}
        return ["va"], ["pf"]
    elseif model_type <: NFAPowerModel
        return [], ["pf"]
    elseif model_type <: DCPLLPowerModel
        return ["va"], ["pf", "pt"]
    elseif model_type <: LPACCPowerModel
        return ["va", "phi"], ["pf", "pt", "qf", "qt", "cs"]
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

"identifyall the variables names"
function variable_names(model_type::DataType)
    if model_type <: Union{DCPPowerModel, DCMPPowerModel}
        return ["va"], ["pf"], ["pg"]
    elseif model_type <: NFAPowerModel
        return [], ["pf"], ["pg"]
    elseif model_type <: DCPLLPowerModel
        return ["va"], ["pf", "pt"], ["pg"]
    elseif model_type <: LPACCPowerModel
        return ["va", "phi"], ["pf", "pt", "qf", "qt", "cs"], ["pg", "qg"]
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