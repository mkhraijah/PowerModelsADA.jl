###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

"""
    solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, kwargs...)

Solve OPF problem using fully distributed algorithm. The distributed algorithm is defined by the build_method and update_method.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follow:
    - initialize_method::Function : initliize distributed algorithm parameters and shared variables
    - update_method::Function : update the algorithm after each iteration
    - build_method::Function : problem formulation
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- verbose::Bool=true : print mismatch after each iteration and result summary 
- print_optimizer_info::Bool=false : print local optimization info from the solver
- kwargs = algorithm parameters

"""
function solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Bool=true, print_optimizer_info::Bool=false, kwargs...)

    ## obtain areas idx
    areas_id = get_areas_id(data)

    ## decompose the system into subsystems
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    ## initilize distributed power model parameters
    for i in areas_id
        dopf_method.initialize_method(data_area[i], model_type; tol=tol, max_iteration=max_iteration, kwargs...)
    end

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        info = @capture_out begin
            ## solve local problem and update solution
            for i in areas_id
                dopf_method.update_method(data_area[i])
                solve_local!(data_area[i], model_type, optimizer, dopf_method.build_method)
                # if model_type <: _PM.AbstractSDPWRMModel
                #     solve_local!(data_area[i], model_type, optimizer, build_method)
                # else
                #     result = _PM.solve_model(data_area[i], model_type, optimizer, build_method)
                #     update_solution!(data_area[i], result["solution"], model_type)
                # end
            end
        end

        if print_optimizer_info
            println(info)
        end

        ## share solution with neighbors, the shared data is first obtained to facilitate distributed implementation  
        for i in areas_id # sender subsystem
            for j in areas_id # receiver subsystem
                if i != j && string(i) in keys(data_area[j]["shared_primal"])
                    shared_data = send_shared_data(i, j, data_area[i], serialize = false)
                    receive_shared_data!(i, shared_data, data_area[j])
                end
            end
        end

        ## calculate mismatches and update convergance flags
        for i in areas_id
            calc_mismatch!(data_area[i], mismatch_method)
            update_flag_convergance!(data_area[i], tol)
            update_iteration!(data_area[i])
        end

        ## check global convergance and update iteration counters
        flag_convergance = update_global_flag_convergance(data_area)
        iteration += 1

        if verbose
            mismatch = calc_global_mismatch(data_area)
            println("Iteration = $iteration, mismatch = $mismatch")
            if flag_convergance
                println("Consistency achived within $tol")
                println("Number of iterations = $iteration")
                objective = calc_dist_gen_cost(data_area)
                println("Objective = $objective")
            end
        end
    end

    return data_area
end

"""
    solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, model_type::DataType; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, kwargs...)

Solve OPF problem using distributed algorithm with central coordinator. The distributed algorithm is defined by the dopf_method Module.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- model_type::DataType : module contains the distributed algorithm methods as follow:
    - initialize_method_coordinator::Function : initliize distributed algorithm parameters and shared variables
    - initialize_method_local::Function : initliize distributed algorithm parameters and shared variables
    - update_method_coordinator::Function : update the algorithm after each iteration
    - update_method_local::Function : update the algorithm after each iteration
    - build_method_coordinator::Function : problem formulation
    - build_method_local::Function : problem formulation
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4, 
- max_iteration::Int64=1000, 
- verbose = true,
- print_optimizer_info::Bool=false : print local optimization info from the solver
- kwargs = distributed algorithm parameters

"""
function solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, 
    dopf_method::Module; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Bool=true, print_optimizer_info::Bool=false, kwargs...)

    ## obtain and arrange areas id
    arrange_areas_id!(data)
    areas_id = get_areas_id(data)


    ## decompose the system into subsystems
    data_coordinator = decompose_coordinator(data)
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end
    
    ## initilize distributed power model parameters
    dopf_method.initialize_method_coordinator(data_coordinator, model_type; tol=tol, max_iteration=max_iteration, kwargs...)
    for i in areas_id
        dopf_method.initialize_method_local(data_area[i], model_type; tol=tol, max_iteration=max_iteration, kwargs...)
    end

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false
        
        # update local and coordinator dual variable
        dopf_method.update_method_coordinator(data_coordinator)
        for i in areas_id
            dopf_method.update_method_local(data_area[i])
        end
        info = @capture_out begin
            # solve local problem and update solution
            for i in areas_id
                solve_local!(data_area[i], model_type, optimizer, dopf_method.build_method_local)
            end
        end

        if print_optimizer_info
            println(info)
        end

        # share solution with the coordinator and add noise to the shared data
        for i in areas_id # sender subsystem
            shared_data = send_shared_data(i, i, data_area[i], serialize = false)
            receive_shared_data!(i, shared_data, data_coordinator)
        end

        info = @capture_out begin
            # solve coordinator problem 
            solve_local!(data_coordinator, model_type, optimizer, dopf_method.build_method_coordinator)
        end

        if print_optimizer_info
            println(info)
        end

        # share coordinator solution with areas and add noise to the shared data
        for i in areas_id # sender subsystem
            shared_data = send_shared_data(0, i, data_coordinator, serialize = false)
            receive_shared_data!(0, shared_data, data_area[i])
        end

        # calculate mismatches and update convergance flags
        for i in areas_id
            calc_mismatch!(data_area[i], mismatch_method)
            update_flag_convergance!(data_area[i], tol)
            update_iteration!(data_area[i])
        end

        calc_mismatch!(data_coordinator, mismatch_method)
        update_flag_convergance!(data_coordinator, tol)
        update_iteration!(data_coordinator)


        # check global convergance and update iteration counters
        flag_convergance = data_coordinator["flag_convergance"]
        iteration += 1

        # print solution
        if verbose
            mismatch = data_coordinator["mismatch"]["0"]
            println("Iteration = $iteration, mismatch = $mismatch")
            if flag_convergance
                println("Consistency achived within $tol")
                println("Number of iterations = $iteration")
                objective = calc_dist_gen_cost(data_area)
                println("Objective = $objective")
            end
        end
    end

    data_area[0] = data_coordinator
    return data_area
end

"initialize dopf parameters"
function initialize_dopf_parameters!(data::Dict{String, <:Any}; tol::Float64=1e-4, max_iteration::Int64=1000)
    data["iteration"] = Int64(0)
    data["flag_convergance"] = false
    data["mismatch"] = Dict{String, Any}()
    data["tol"] = tol
    data["max_iteration"] = max_iteration
end

"update the area data and local shared variables after obtaining a solution at each iteraton"
function update_solution!(data::Dict{String, <:Any}, solution::Dict{String, <:Any}, model_type::DataType)
    _PM.update_data!(data, solution)
    update_shared_primal!(data, model_type)
end

"update primal variables after obtaining a solution at each iteraton"
function update_shared_primal!(data::Dict{String,<:Any}, model_type::DataType)
    area_id = string(get_area_id(data))
    bus_variables_name, branch_variables_name =  variable_shared_names(model_type)
    shared_primal = data["shared_primal"][area_id]

    for var in bus_variables_name
        for idx in keys(shared_primal[var])
            shared_primal[var][idx] = data["bus"][idx][var]
        end
    end
    for var in branch_variables_name
        for idx in keys(shared_primal[var])
            shared_primal[var][idx] = data["branch"][idx][var]
        end
    end
end

function solve_local!(data::Dict{String, <:Any}, model_type::DataType, optimizer, build_method::Function)
    pm = _PM.instantiate_model(data, model_type, build_method)
    result = _PM.optimize_model!(pm, relax_integrality=false, optimizer=optimizer, solution_processors=[])
    _PM.update_data!(data, result["solution"])

    area_id = string(get_area_id(data))
    shared_primal = data["shared_primal"][area_id]

    for var in keys(shared_primal)
        for idx in keys(shared_primal[var])
            shared_primal[var][idx] =  JuMP.value(_var(pm, var, idx))
        end
    end
end

"""
    calc_mismatch!(data::Dict{String, <:Any},method::String="norm"; p::Int64=2)
calculate the mismatch using p-norm and return the area data dictionary with the mismatch as seen by the area.
"""
function calc_mismatch!(data::Dict{String, <:Any}, method::String="norm"; p::Int64=2 )
    area_id = string(get_area_id(data)) 
    primal_variable = data["shared_primal"]

    mismatch = Dict{String, Any}([
        area => Dict{String, Any}([
            variable => Dict{String, Any}([
                idx => primal_variable[area_id][variable][idx] - primal_variable[area][variable][idx]
            for idx in keys(primal_variable[area][variable])]) 
        for variable in keys(primal_variable[area])]) 
    for area in keys(primal_variable) if area != area_id ])

    if method == "norm"
        mismatch[area_id] = LinearAlgebra.norm([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]], p)
    elseif method == "max" || method == "maximum"
        mismatch[area_id] = LinearAlgebra.maximum([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]])
    end

    data["mismatch"] = mismatch
end

"calculate the global mismatch based on local mismatch"
function calc_global_mismatch(data_area::Dict{Int, <:Any}, method::String="norm"; p::Int64=2) 
    if method == "norm"
        LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area)], p)
    elseif method == "max" || method == "maximum"
        LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area)])
    end
end

"check the shared variables of a local area are within tol"
function update_flag_convergance!(data::Dict{String, <:Any}, tol::Float64)
    area_id = string(data["area"])
    mismatch = data["mismatch"][area_id]
    data["flag_convergance"] = mismatch < tol
end

"check the flag convergence for all areas and return a global variables"
function update_global_flag_convergance(data_area::Dict{Int, <:Any})
    flag_convergance = reduce( & , [data_area[i]["flag_convergance"] for i in keys(data_area)])
    return flag_convergance
end

"update iteration"
function update_iteration!(data::Dict{String, <:Any})
    data["iteration"] += 1
end

"""
arrange area id from 1 to number of areas
this step is necessary when having area number 0 and using central coordinator
"""
function arrange_areas_id!(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    new_areas_id = collect(1:length(areas_id))
    area_id_lookup = Dict(areas_id[i] => i for i in new_areas_id)
    for (i,bus) in data["bus"]
        bus["area"] = area_id_lookup[bus["area"]]
    end
end

"return JuMP variable object from PowerModel object"
function _var(pm::AbstractPowerModel, key::String, idx::String)
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

"helper functions to handle area ids, local buses, neighbor buses"
get_areas_id(data::Dict{String, <:Any}) = unique([bus["area"] for (i, bus) in data["bus"]])

get_areas_id(pm::AbstractPowerModel) = get_areas_id(pm.data)

get_area_id(data::Dict{String, <:Any}) = get(data,"area", NaN)

get_area_id(pm::AbstractPowerModel) = get_area_id(pm.data)

get_local_bus(data::Dict{String, <:Any}, area::Int) = Vector{Int64}([bus["bus_i"] for (i,bus) in data["bus"] if bus["area"] == area])

get_local_bus(pm::AbstractPowerModel, area::Int) = get_local_bus(pm.data, area)

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

get_neighbor_bus(pm::AbstractPowerModel, local_bus::Vector) = get_neighbor_bus(pm.data, local_bus)

get_neighbor_bus(data::Dict{String, <:Any}, area::Int) = get_neighbor_bus(data, get_local_bus(data,area))

get_neighbor_bus(pm::AbstractPowerModel, area::Int) = get_neighbor_bus(pm.data, area)

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

get_shared_component(pm::AbstractPowerModel, area::Int64) = get_shared_component(pm.data, area)

function get_shared_component(data::Dict{String, <:Any})
    area = get_area_id(data)
    get_shared_component(data, area)
end

get_shared_component(pm::AbstractPowerModel) = get_shared_component(pm.data)