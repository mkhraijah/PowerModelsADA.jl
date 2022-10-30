###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

"""
    solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, print_level::Int64=1, kwargs...)

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
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- kwargs = algorithm parameters
"""
function solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)
    # get areas ids
    areas_id = get_areas_id(data)

    # decompose the system into subsystems
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    # initilize distributed power model parameters
    for i in areas_id
        dopf_method.initialize_method(data_area[i], model_type; kwargs...)
    end

    # get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    # initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    # start iteration
    while iteration < max_iteration && !flag_convergance

        # solve local problem and update solution
        info = @capture_out begin
            Threads.@threads for i in areas_id
                result = solve_model(data_area[i], model_type, optimizer, dopf_method.build_method, solution_processors=dopf_method.post_processors)
                update_data!(data_area[i], result["solution"])
            end
        end

        # share solution with neighbors, the shared data is first obtained to facilitate distributed implementation
        for i in areas_id # sender subsystem
            for j in areas_id # receiver subsystem
                if i != j
                    shared_data = prepare_shared_data(data_area[i], j)
                    receive_shared_data!(data_area[j], shared_data, i)
                end
            end
        end

        # calculate mismatches and update convergance flags
        Threads.@threads for i in areas_id
            dopf_method.update_method(data_area[i])
        end

        # check global convergance and update iteration counters
        flag_convergance = update_global_flag_convergance(data_area)
        iteration += 1

        print_iteration(data_area, print_level, [info])
        
    end

    print_convergance(data_area, print_level)

    return data_area
end


function solve_dopf(data::String, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)
    data = parse_file(data)
    solve_dopf(data, model_type, optimizer, dopf_method; print_level=print_level, kwargs...)
end

"""
    solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, model_type::DataType; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, print_level::Int64=1, kwargs...)

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
- print_level::Int64=1 : print mismatch after each iteration and result summary 
- kwargs = distributed algorithm parameters

"""
function solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)

    ## arrange and get areas id
    arrange_areas_id!(data)
    areas_id = get_areas_id(data)

    ## decompose the system into subsystems
    data_coordinator = decompose_coordinator(data)
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end
    
    ## initilize distributed power model parameters
    dopf_method.initialize_method_coordinator(data_coordinator, model_type; kwargs...)
    for i in areas_id
        dopf_method.initialize_method_local(data_area[i], model_type; kwargs...)
    end

    ## get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && !flag_convergance

        # solve local area problems in parallel
        info1 = @capture_out begin
            Threads.@threads for i in areas_id
                result = solve_model(data_area[i], model_type, optimizer, dopf_method.build_method_local, solution_processors=dopf_method.post_processors_local)
                update_data!(data_area[i], result["solution"])
            end
        end

        # share solution of local areas with the coordinator
        for i in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_area[i], 0, serialize = false)
            receive_shared_data!(data_coordinator, shared_data, i)
        end

        # solve coordinator problem 
        info2 = @capture_out begin
            result = solve_model(data_coordinator, model_type, optimizer, dopf_method.build_method_coordinator, solution_processors=dopf_method.post_processors_coordinator)
            update_data!(data_coordinator, result["solution"])
        end

        # share coordinator solution with local areas
        for i in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_coordinator, i, serialize = false)
            receive_shared_data!(data_area[i], shared_data, 0)
        end

        # update local areas and coordinator problems after
        dopf_method.update_method_coordinator(data_coordinator)
        for i in areas_id
            dopf_method.update_method_local(data_area[i])
        end

        # check global convergance and update iteration counters
        flag_convergance = data_coordinator["counter"]["flag_convergance"]
        iteration += 1

        # print solution
        print_iteration(data_area, print_level, [info1; info2])

    end

    data_area[0] = data_coordinator
    print_convergance(data_area, print_level)

    return data_area
end

function solve_dopf_coordinated(data::String, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)
    data = parse_file(data)
    solve_dopf_coordinated(data, model_type, optimizer, dopf_method; print_level=print_level, kwargs...)
end

"initialize dopf parameters"
function initialize_dopf!(data::Dict{String, <:Any}, model_type::DataType; kwargs...)
    # options
    data["option"] = Dict{String, Any}()
    data["option"]["tol"] = get(kwargs, :tol, 1e-4)
    data["option"]["max_iteration"] = get(kwargs, :max_iteration, 1000)
    data["option"]["mismatch_method"] = get(kwargs, :mismatch_method, "norm")
    data["option"]["model_type"] = model_type

    # counters
    data["counter"] = Dict{String, Any}()
    data["counter"]["iteration"] = Int64(0)
    data["counter"]["flag_convergance"] = false

    # mismatch 
    data["mismatch"] = Dict{String, Any}()

    # last solution 
    data["solution"] = initialize_all_variable(data, model_type)

    # previous solutions 
    save_data = get(kwargs, :save_data, ["solution"])
    if !isempty(save_data)
        data["previous_solution"] = Dict{String, Any}([str=>Vector{Dict}() for str in save_data])
    end
end

"update the area data solution dictionary"
function update_solution!(pm::AbstractPowerModel, solution::Dict{String, <:Any})
    solution["solution"] = pm.data["solution"]
    for variable in keys(solution["solution"])
        for idx in keys(solution["solution"][variable])
            solution["solution"][variable][idx] = JuMP.value(PowerModelsADA._var(pm, variable, idx))
        end
    end
end

"update primal variables after obtaining a solution at each iteraton"
function update_shared_variable!(pm::AbstractPowerModel, solution::Dict{String, <:Any})
    solution["shared_variable"] = pm.data["shared_variable"]
    for area in keys(solution["shared_variable"])
        for var in keys(solution["shared_variable"][area])
            for idx in keys(solution["shared_variable"][area][var])
                solution["shared_variable"][area][var][idx] = solution["solution"][var][idx]
            end
        end
    end
end

"save last solution in previous_solutions vector"
function save_solution!(data::Dict{String, <:Any})
    if haskey(data, "previous_solution")
        for str in keys(data["previous_solution"])
            push!(data["previous_solution"][str], deepcopy(data[str]))
        end
    end
end

"update iteration"
function update_iteration!(data::Dict{String, <:Any})
    data["counter"]["iteration"] += 1
end

"""
    calc_mismatch!(data::Dict{String, <:Any},method::String="norm"; p::Int64=2)
calculate the mismatch using p-norm and return the area data dictionary with the mismatch as seen by the area.
"""
function calc_mismatch!(data::Dict{String, <:Any}; p::Int64=2 )
    area_id = string(get_area_id(data)) 
    mismatch_method = data["option"]["mismatch_method"]
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]

    mismatch = Dict{String, Any}([
        area => Dict{String, Any}([
            variable => Dict{String, Any}([
                idx => shared_variable_local[area][variable][idx] - shared_variable_received[area][variable][idx]
            for idx in keys(shared_variable_local[area][variable])]) 
        for variable in keys(shared_variable_local[area])]) 
    for area in keys(shared_variable_local) if area != area_id && area in keys(shared_variable_received) ])

    if mismatch_method == "norm"
        mismatch[area_id] = LinearAlgebra.norm([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]], p)
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        mismatch[area_id] = LinearAlgebra.maximum([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]])
    end

    data["mismatch"] = mismatch
end

"check the shared variables of a local area are within tol"
function update_flag_convergance!(data::Dict{String, <:Any})
    area_id = string(data["area"])
    tol = data["option"]["tol"]
    mismatch = data["mismatch"][area_id]
    data["counter"]["flag_convergance"] = mismatch < tol
end

"calculate the global mismatch based on local mismatch"
function calc_global_mismatch(data_area::Dict{Int, <:Any}; p::Int64=2) 
    mismatch_method = first(data_area)[2]["option"]["mismatch_method"]
    if mismatch_method == "norm"
        LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0], p)
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
    end
end

"check the flag convergence for all areas and return a global variables"
function update_global_flag_convergance(data_area::Dict{Int, <:Any}, central::Bool=true)
    if central
        mismatch = calc_global_mismatch(data_area)
        tol =first(data_area)[2]["option"]["tol"]
        return mismatch <= tol
    else
        return reduce( &, [data_area[i]["counter"]["flag_convergance"] for i in keys(data_area)])
    end
end 

function print_iteration(data::Dict, print_level::Int64, info_list::Vector=[])
    if print_level > 0
        iteration = first(data)[2]["counter"]["iteration"]
        mismatch = calc_global_mismatch(data)
        println("Iteration = $iteration, mismatch = $mismatch")

        if print_level > 1
            for info in info_list
                println(info)
            end
        end
    end
end

function print_convergance(data::Dict, print_level::Int64)
    if print_level > 0
        iteration = first(data)[2]["counter"]["iteration"]
        mismatch = calc_global_mismatch(data)
        tol =first(data)[2]["option"]["tol"]
        flag_convergance = mismatch <= tol
        if flag_convergance
            println("*******************************************************")
            println("")
            println("Consistency achived within $tol mismatch tolerence")
            println("Number of iterations = $iteration")
            
            objective = calc_dist_gen_cost(data)
            println("Objective function value = $objective")
            println("")
            println("*******************************************************")
        else
            println("*******************************************************")
            println("")
            println("Consistency did not achived within $tol mismatch tolerence and $iteration iteration")
            println("Shared variables mismatch = $mismatch")
            println("")
            println("*******************************************************")
        end
    end
end