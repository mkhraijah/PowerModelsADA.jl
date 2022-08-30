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


"helper functions to handle area ids, local buses, neighbor buses"
get_areas_id(data::Dict{String, <:Any}) = unique([bus["area"] for (i, bus) in data["bus"]])

get_areas_id(pm::AbstractPowerModel) = unique([bus["area"] for (i, bus) in pm.data["bus"]])

get_area_id(data::Dict{String, <:Any}) = get(data,"area", NaN)

get_area_id(pm::AbstractPowerModel) = get(pm.data,"area", NaN)

get_local_bus(data::Dict{String, <:Any},area::Int) = Vector{Int64}([bus["bus_i"] for (i,bus) in data["bus"] if bus["area"] == area])

get_local_bus(pm::AbstractPowerModel,area::Int) = Vector{Int64}([bus["bus_i"] for (i,bus) in pm.data["bus"] if bus["area"] == area])

function get_neighbor_bus(data::Dict{String, <:Any}, local_bus::Vector)
    neighbor_bus = Vector{Int64}()
    for (i,branch) in data["branch"]
        if branch["f_bus"] in local_bus && !(branch["t_bus"] in local_bus)
            push!(neighbor_bus,branch["t_bus"])
        elseif !(branch["f_bus"] in local_bus) && branch["t_bus"] in local_bus
            push!(neighbor_bus,branch["f_bus"])
        end
    end
    return neighbor_bus
end

get_neighbor_bus(pm::AbstractPowerModel, area::Int) = get_neighbor_bus(pm.data, area)

get_neighbor_bus(data::Dict{String, <:Any}, area::Int) = get_neighbor_bus(data, get_local_bus(data,area))

function get_areas_bus(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    areas_bus = Dict{Int64, Vector{Int64}}()
    for i in areas_id
        areas_bus[i] = [bus["bus_i"] for (j,bus) in data["bus"] if bus["area"]==i]
    end
    return areas_bus
end

get_areas_bus(pm::AbstractPowerModel) = get_areas_bus(pm.data)

"get the shared buses and branches between defined area and all other areas"
function get_shared_component(data::Dict{String, <:Any}, area::Int64)
    areas_id = get_areas_id(data)
    areas_bus = get_areas_bus(data)
    shared_branch = Dict{Int64, Any}()
    shared_bus = Dict{Int64, Any}()
    for i in areas_id
        if i != area
            shared_branch[i] = unique([parse(Int64,j) for (j,branch) in data["branch"] if (branch["f_bus"] in areas_bus[i] && branch["t_bus"] in areas_bus[area]) || (branch["f_bus"] in areas_bus[area] && branch["t_bus"] in areas_bus[i]) ])
        else
            shared_branch[i] = unique([parse(Int64,j) for (j,branch) in data["branch"] if xor(branch["f_bus"] in areas_bus[i], branch["t_bus"] in areas_bus[i]) ])
        end
            shared_bus[i] = unique(vcat([branch["f_bus"] for (j,branch) in data["branch"] if parse(Int64,j) in shared_branch[i]], [branch["t_bus"] for (j,branch) in data["branch"] if parse(Int64,j) in shared_branch[i]] ))
    end
    return shared_bus, shared_branch
end

function get_shared_component(data::Dict{String, <:Any})
    area = get_area_id(data)
    get_shared_component(data, area)
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
    Central_solution = _PM.run_opf(data, model_type, optimizer)

    # Calculate objective function
    Obj_dist = calc_dist_gen_cost(data_area::Dict{Int, <:Any})
    Obj_cent = Central_solution["objective"]

    # Calculate optimility gap
    Relative_Error = (Obj_dist - Obj_cent)/ Obj_cent * 100
    return Relative_Error
end


function _var(pm::AbstractPowerModel, key::String, idx::String )
    bus_variables_name, branch_variables_name = variable_shared_names(typeof(pm))
    idx = parse(Int64,idx)
    if key in bus_variables_name
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