###############################################################################
#               Helper methods for all distributed algorithms                 #
###############################################################################

"partition a system into n areas"
function partition_system!(data::Dict, n::Int64; configuration::Symbol=:edge_cut, print_info::Bool=false)
 
    nbus = length(data["bus"])
    nbranch = length(data["branch"])
    bus_index = [x.second["index"] for x in data["bus"]]
    branch_index = [x.second["index"] for x in data["branch"]]

    sort!(bus_index)
    sort!(branch_index)

    W = zeros(nbus,nbranch)

    for (i,branch) in data["branch"]
        f_bus = findfirst(x->x==branch["f_bus"], bus_index)
        t_bus = findfirst(x->x==branch["t_bus"], bus_index)
        indx = branch["index"]
        W[f_bus,indx] = 1
        W[t_bus,indx] = 1
    end
    W = sparse(W)
    h = KaHyPar.HyperGraph(W)

    info = @capture_out begin
        partitions = KaHyPar.partition(h, n, configuration=configuration)
    end
    partitions = Dict([bus_index[i]=>partitions[i] for i in 1:nbus])

    for (i,bus) in data["bus"]
        bus["area"] = partitions[bus["index"]]
    end

    arrange_areas_id!(data)
    if print_info
        println(info)
    end
end


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
    Obj_dist = calc_dist_gen_cost(data_area::Dict{Int, <:Any})
    Obj_cent = Central_solution["objective"]
    # Calculate optimility gap
    Relative_Error = (Obj_dist - Obj_cent)/ Obj_cent * 100
    return Relative_Error
end
