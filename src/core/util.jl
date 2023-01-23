###############################################################################
#               Helper methods for all distributed algorithms                 #
###############################################################################

# partition_system! function is removed due to a compiling issue of the KaHyPar package with Windows 
# """
#     partition_system!(data::Dict, n::Int64; configuration::Symbol=:edge_cut, print_info::Bool=false)

# Partition a system into n areas using KaHyPar partition algorithm

# # Arguments:
# - data::Dict{String, <:Any} : dictionary contains case in PowerModel format
# - n::Int : number of areas
# - configuration::Symbol=:edge_cut : partition meteric (:edge_cut or :connectivity)
# - print_info::Bool=false : print partition algorithm information
# """
# function partition_system!(data::Dict, n::Int64; configuration::Symbol=:edge_cut, print_info::Bool=false)
#     nbus = length(data["bus"])
#     nbranch = length(data["branch"])
#     bus_index = [x.second["index"] for x in data["bus"]]
#     branch_index = [x.second["index"] for x in data["branch"]]

#     sort!(bus_index)
#     sort!(branch_index)

#     W = zeros(nbus,nbranch)

#     for (i,branch) in data["branch"]
#         f_bus = findfirst(x->x==branch["f_bus"], bus_index)
#         t_bus = findfirst(x->x==branch["t_bus"], bus_index)
#         indx = branch["index"]
#         W[f_bus,indx] = 1
#         W[t_bus,indx] = 1
#     end
#     W = SparseArrays.sparse(W)
#     h = KaHyPar.HyperGraph(W)

#     info = @capture_out begin
#         partitions = KaHyPar.partition(h, n, configuration=configuration)
#     end
#     partitions = Dict([bus_index[i]=>partitions[i]+1 for i in 1:nbus])

#     for (i,bus) in data["bus"]
#         bus["area"] = partitions[bus["index"]]
#     end

#     if print_info
#         println(info)
#     end
# end

"calculate distributed solution operation cost"
function calc_dist_gen_cost(data_area::Dict{Int, <:Any})
    gen_cost = 0
    # Calculate objective function
    for i in keys(data_area)
        gen_cost += _PM.calc_gen_cost(data_area[i])
    end
    return gen_cost
end

"compare the distributed algorithm solution with PowerModels centralized solution"
function compare_solution(data::Dict{String, <:Any}, data_area::Dict{Int, <:Any}, model_type::DataType, optimizer)
    # Solve Centralized OPF
    Central_solution = _PM.solve_opf(data, model_type, optimizer)
    # Calculate objective function
    Obj_distributed = calc_dist_gen_cost(data_area::Dict{Int, <:Any})
    Obj_centeral = Central_solution["objective"]
    # Calculate optimality gap
    Relative_Error = abs(Obj_distributed - Obj_centeral)/ Obj_centeral * 100
    return Relative_Error
end


"get the number of variables in each area"
function calc_number_areas_variables(data::Dict{String, <:Any}, model_type::DataType)
    areas_id = get_areas_id(data)
    data_area = Dict{Int64, Any}()
    num_variables = Dict{Int64, Any}()
    for area in areas_id
        data_area = decompose_system(data, area)
        num_variables[area] = calc_number_all_variables(data_area, model_type)
    end
    return num_variables
end

"get the number of shared variable in a area"
function calc_number_shared_variables(data::Dict{String, <:Any}, model_type::DataType)
    areas_id = get_areas_id(data)
    area_id = get_area_id(data)
    shared_variable = initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", "flat")
    num_variables = calc_number_variables(shared_variable)
    return num_variables
end

"get the number of shared variable in all area"
function calc_number_system_shared_variables(data::Dict{String, <:Any}, model_type::DataType)
    areas_id = get_areas_id(data)
    num_variables = Dict()
    for area in areas_id
        data_area = decompose_system(data, area)
        num_variables[area] = calc_number_shared_variables(data_area, model_type)
    end
    return num_variables
end

"get the number of variables in an area"
function calc_number_all_variables(data::Dict{String, <:Any}, model_type::DataType)
    variables = initialize_all_variable(data, model_type)
    num_variables = calc_number_variables(variables)
    return num_variables
end

"get the number of variables"
function calc_number_variables(data::Dict{String, <:Any})
    num_variables = 0 
    for (key,val) in data
        if isa(val, Dict{String, <:Any})
            num_variables += calc_number_variables(val)
        else
            num_variables += length(val)
        end
    end
    return num_variables
end