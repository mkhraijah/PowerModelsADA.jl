###############################################################################
#               Helper methods for all distributed algorithms                 #
###############################################################################

"partition a system into p subsystem using spactural clustering"
function partition_system!(data::Dict, p::Int64; init=[], correct_neg_susceptance::Bool=false)

    neg_b_branch = findall(x->x["br_x"]<0, data["branch"])
    if !isempty(neg_b_branch)
        if !correct_neg_susceptance
            error("Branch/s $neg_b_branch has/have negative susceptance, try to enforce postive susceptance by setting correct_neg_susceptance = true")
        end
    end

    nbus = length(data["bus"])

    bus_index = [x.second["index"] for x in data["bus"]]
    sort!(bus_index)

    W = zeros(nbus,nbus)
    D = zeros(nbus,nbus)

    for (i,branch) in data["branch"]
        f_bus = findfirst(x->x==branch["f_bus"], bus_index)
        t_bus = findfirst(x->x==branch["t_bus"], bus_index)
        br_x = abs(branch["br_x"])
        W[f_bus,t_bus] = 1/br_x
        W[t_bus,f_bus] = 1/br_x
    end

    for i in 1:nbus
        D[i,i] = sum(W[i,:])
    end

    L = LinearAlgebra.I - sqrt.(inv(D)) * W * sqrt.(inv(D))
    S,U = LinearAlgebra.eigen(L)
    if isempty(init)
        C = kmeans(U[:,1:p]',p).assignments
    else
        order = Vector{Int64}()
        for i in init
            append!(order, findfirst(isequal(i), bus_index))
        end
        C = kmeans(U[:,1:p]',p, init=order).assignments
    end

    for i in 1:nbus
        data["bus"]["$(bus_index[i])"]["area"] = C[i]
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
function compare_solution(data::Dict{String, <:Any}, data_area::Dict{Int, <:Any}, model_type, optimizer)
    # Solve Centralized OPF
    Central_solution = _PM.solve_opf(data, model_type, optimizer)
    # Calculate objective function
    Obj_dist = calc_dist_gen_cost(data_area::Dict{Int, <:Any})
    Obj_cent = Central_solution["objective"]
    # Calculate optimility gap
    Relative_Error = (Obj_dist - Obj_cent)/ Obj_cent * 100
    return Relative_Error
end