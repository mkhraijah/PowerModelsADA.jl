###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

"""
    solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Int64=1, kwargs...)

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
- verbose::Int64=1 : print mismatch after each iteration and result summary 
- kwargs = algorithm parameters
"""
# mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000,
function solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; verbose::Int64=1, kwargs...)

    # obtain areas ids
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

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        # solve local problem and update solution
        info = @capture_out begin
            Threads.@threads for i in areas_id
                result = _PM.solve_model(data_area[i], model_type, optimizer, dopf_method.build_method, solution_processors=[update_solution!, update_shared_variable!] )
                update_data!(data_area[i], result["solution"])
            end
        end

        # share solution with neighbors, the shared data is first obtained to facilitate distributed implementation
        for i in areas_id # sender subsystem
            for j in areas_id # receiver subsystem
                if i != j
                    shared_data = prepare_shared_data(data_area[i], j, serialize = false)
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

        print_iteration(data_area, verbose, [info])
        
    end

    return data_area
end


function solve_dopf(data::String, model_type::DataType, optimizer, dopf_method::Module; verbose::Int64=1, kwargs...)

    data = parse_file(data)
    solve_dopf(data, model_type, optimizer, dopf_method; verbose=verbose, kwargs...)

end

"""
    solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, model_type::DataType; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Int64=1, kwargs...)

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
- verbose::Int64=1 : print mismatch after each iteration and result summary 
- kwargs = distributed algorithm parameters

"""
function solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, 
    dopf_method::Module; verbose::Int64=1, kwargs...)

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
    while iteration < max_iteration && flag_convergance == false

        # solve local area problems in parallel
        info1 = @capture_out begin
            Threads.@threads for i in areas_id
                result = _PM.solve_model(data_area[i], model_type, optimizer, dopf_method.build_method_local, solution_processors=[update_solution!, update_shared_variable!] )
                update_data!(data_area[i], result["solution"])
            end
        end

        # share solution of local areas with the coordinator
        for i in areas_id # sender subsystem
            shared_data = PowerModelsADA.prepare_shared_data(data_area[i], 0, serialize = false)
            PowerModelsADA.receive_shared_data!(data_coordinator, shared_data, i)
        end

        # solve coordinator problem 
        info2 = @capture_out begin
            result = _PM.solve_model(data_coordinator, model_type, optimizer, dopf_method.build_method_coordinator, solution_processors=[update_solution!, update_shared_variable!] )
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
        print_iteration(data_area, verbose, [info1; info2])

    end

    data_area[0] = data_coordinator
    return data_area
end

function solve_dopf_coordinated(data::String, model_type::DataType, optimizer, dopf_method::Module; verbose::Int64=1, kwargs...)

    data = parse_file(data)
    solve_dopf_coordinated(data, model_type, optimizer, dopf_method; verbose=verbose, kwargs...)
end

"initialize dopf parameters"
function initialize_dopf_parameters!(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    # options
    data["option"] = Dict{String, Any}()
    data["option"]["tol"] = get(kwargs, :tol, 1e-4)
    data["option"]["max_iteration"] =  get(kwargs, :max_iteration, 1000)
    data["option"]["mismatch_method"] =  get(kwargs, :mismatch_method, "norm")
    data["option"]["model_type"] =  model_type

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
    if haskey(pm.data, "solution")
        solution["solution"] = pm.data["solution"]
        for variable in keys(solution["solution"])
            for idx in keys(solution["solution"][variable])
                solution["solution"][variable][idx] = JuMP.value(PowerModelsADA._var(pm, variable, idx))
            end
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

# "solve local problem and update the solution and shared variables dictionaries"
# function solve_local!(data::Dict{String, <:Any}, model_type::DataType, optimizer, build_method::Function)

#     pm = _PM.instantiate_model(data, model_type, build_method)
#     result = _PM.optimize_model!(pm, relax_integrality=false, optimizer=optimizer, solution_processors=[])
#     _PM.update_data!(data, result["solution"])

#     if haskey(data, "solution")
#         solution = data["solution"]
#         for variable in keys(solution)
#             for idx in keys(solution[variable])
#                 solution[variable][idx] = JuMP.value(_var(pm, variable, idx))
#             end
#         end
#     end

#     area_id = string(get_area_id(data))
#     shared_variable = data["shared_variable"][area_id]
#     for variable in keys(shared_variable)
#         for idx in keys(shared_variable[variable])
#             shared_variable[variable][idx] = JuMP.value(_var(pm, variable, idx))
#         end
#     end

# end

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
    shared_variable_received = data["received_shared_variable"]

    mismatch = Dict{String, Any}([
        area => Dict{String, Any}([
            variable => Dict{String, Any}([
                idx => shared_variable_local[area][variable][idx] - shared_variable_received[area][variable][idx]
            for idx in keys(shared_variable_local[area][variable])]) 
        for variable in keys(shared_variable_local[area])]) 
    for area in keys(shared_variable_local) if area != area_id ])

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
        LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area)], p)
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area)])
    end
end

"check the flag convergence for all areas and return a global variables"
update_global_flag_convergance(data_area::Dict{Int, <:Any}) = reduce( & , [data_area[i]["counter"]["flag_convergance"] for i in keys(data_area)])

function print_iteration(data::Dict, verbose::Int64, info_list::Vector=[])

    if verbose > 0
        iteration = first(data)[2]["counter"]["iteration"]
        mismatch = calc_global_mismatch(data)
        flag_convergance = mismatch <= first(data)[2]["option"]["tol"]
        println("Iteration = $iteration, mismatch = $mismatch")
        if flag_convergance
            tol =first(data)[2]["option"]["tol"]
            println("*******************************************************")
            println("")
            println("Consistency achived within $tol mismatch tolerence")
            println("Number of iterations = $iteration")
            objective = calc_dist_gen_cost(data)
            println("Objective function value = $objective")
            println("")
            println("*******************************************************")
        end
    end

    if verbose > 1
        for info in info_list
            println(info)
        end
    end
end