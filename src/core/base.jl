###############################################################################
#              Base method for all distributed OPF algorithms                 #
###############################################################################


""
function solve_pmada_model(data::Dict{String,<:Any}, model_type::Type, optimizer, build_method;
        ref_extensions=[], solution_processors=[], relax_integrality=false,
        multinetwork=false, multiconductor=false, kwargs...)

    if multinetwork != _IM.ismultinetwork(data)
        model_requirement = multinetwork ? "multi-network" : "single-network"
        data_type = _IM.ismultinetwork(data) ? "multi-network" : "single-network"
    end

    if multiconductor != ismulticonductor(data)
        model_requirement = multiconductor ? "multi-conductor" : "single-conductor"
        data_type = ismulticonductor(data) ? "multi-conductor" : "single-conductor"
    end

    pm = instantiate_pmada_model(data, model_type, build_method; ref_extensions=ref_extensions, kwargs...)

    result = optimize_model!(pm, relax_integrality=relax_integrality, optimizer=optimizer, solution_processors=solution_processors)

    return result
end

""
function instantiate_pmada_model(data::Dict{String,<:Any}, model_type::Type, build_method; kwargs...)
    return _IM.instantiate_model(data, model_type, build_method, ref_add_core!, _pmada_global_keys, pm_it_sym; kwargs...)
end

""
function build_pmada_ref(data::Dict{String,<:Any}; ref_extensions=[])
    return _IM.build_ref(data, ref_add_core!, _pmada_global_keys, pm_it_name; ref_extensions=ref_extensions)
end

"""
    solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using fully distributed algorithm.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize_method::Function : initialize the algorithm parameters and shared variables
    - update_method::Function : update the algorithm after each iteration
    - build_method::Function : problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- multiprocessors::Bool=false : enable multiprocessors using available workers. Multiprocessors feature requires loading the PowerModelsADA and the optimizer packages on all the processors using @everywhere using <package_name>.
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    # arrange and get areas id
    arrange_areas_id!(data)
    areas_id = get_areas_id(data)
    diameter = get_diameter(data)

    if length(areas_id) < 2
        error("Number of areas is less than 2, at least 2 areas is needed")
    end

    # decompose the system into subsystems
    data_area = Dict{Int64, Any}()
    for area in areas_id
        data_area[area] = decompose_system(data, area)
    end

    solve_dopf(data_area, model_type, optimizer, dopf_method; print_level, multiprocessors=multiprocessors, diameter=diameter, all_areas=areas_id, kwargs...)
end

function solve_dopf(data::String, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    data = parse_file(data)
    solve_dopf(data, model_type, optimizer, dopf_method; print_level=print_level, multiprocessors=multiprocessors, kwargs...)
end

function solve_dopf(data::Dict{Int64, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    if !multiprocessors
        solve_dopf_sp(data, model_type, optimizer, dopf_method; print_level, kwargs...)
    else
        solve_dopf_mp(data, model_type, optimizer, dopf_method; print_level, kwargs...)
    end
end

"""
    solve_dopf_sp(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using fully distributed algorithm on single-processor.
# Arguments:
- data::Dict{Int64, <:Any} : dictionary contains area data in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize_method::Function : initialize the algorithm parameters and shared variables
    - update_method::Function : update the algorithm after each iteration
    - build_method::Function : problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf_sp(data_area::Dict{Int64, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)
    # get areas ids
    areas_id = get_areas_id(data_area)

    # initilize distributed power model parameters
    for area in areas_id
        dopf_method.initialize_method(data_area[area], model_type; kwargs...)
    end

    # get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    # initialize the algorithms global counters
    iteration = 1
    flag_convergence = false

    # start iteration
    while iteration <= max_iteration && !flag_convergence

        # solve local problem and update solution
        info = @capture_out begin
            Threads.@threads for area in areas_id
                result = solve_pmada_model(data_area[area], model_type, optimizer, dopf_method.build_method, solution_processors=dopf_method.post_processors)
                update_data!(data_area[area], result["solution"])
            end
        end

        # share solution with neighbors, the shared data is first obtained to facilitate distributed implementation
        for area in areas_id # sender subsystem
            for neighbor in data_area[area]["neighbors"] # receiver subsystem
                shared_data = prepare_shared_data(data_area[area], neighbor)
                receive_shared_data!(data_area[neighbor], deepcopy(shared_data), area)
            end
        end

        # calculate mismatches and update convergence flags
        Threads.@threads for area in areas_id
            dopf_method.update_method(data_area[area])
        end

        # print solution
        print_iteration(data_area, print_level, [info])

        # check global convergence and update iteration counters
        flag_convergence = update_global_flag_convergence(data_area)
        iteration += 1

    end

    print_convergence(data_area, print_level)
    return data_area
end

"""
    solve_dopf_mp(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using fully distributed algorithm on multiprocessors. Multiprocessors feature requires loading the PowerModelsADA and the optimizer packages on all the processors using @everywhere using <package_name>.
# Arguments:
- data::Dict{Int64, <:Any} : dictionary contains area data in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize_method::Function : initialize the algorithm parameters and shared variables
    - update_method::Function : update the algorithm after each iteration
    - build_method::Function : problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf_mp(data_area::Dict{Int64, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)

    # lookup dictionaries for worker-area pairs
    areas_id = get_areas_id(data_area)
    worker_id = Distributed.workers()
    number_workers = length(worker_id)

    k = 1
    area_worker = Dict()
    for i in areas_id
        if k > number_workers
            k = 1
        end
        area_worker[i] = worker_id[k]
        k += 1 
    end
    worker_area = Dict([i => findall(x -> x==i, area_worker) for i in worker_id if i in values(area_worker)])

    # initiate communication channels 
    comms = Dict(0 => Dict(area => Distributed.RemoteChannel(1) for area in areas_id))
    for area1 in areas_id
        comms[area1] = Dict()
        for area2 in [0; areas_id]
            if area1 != area2 
                comms[area1][area2] = Distributed.RemoteChannel(area_worker[area1])
            end
        end
    end

    # initilize distributed power model parameters
    for area in areas_id
        dopf_method.initialize_method(data_area[area], model_type; kwargs...)
        put!(comms[0][area], data_area[area])
    end

    # get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    # initialize the algorithms global counters
    iteration = 1
    global_flag_convergence = false
    global_counters = Dict{Int64, Any}()

    # share global variables
    Distributed.@everywhere keys(worker_area) begin
        comms = $comms
        areas_id = $areas_id
        worker_area = $worker_area
        area_worker = $area_worker
        dopf_method = $dopf_method
        model_type = $model_type
        optimizer = $optimizer
        area_id = worker_area[myid()]
        data_local = Dict{Int64, Any}(area => take!(comms[0][area]) for area in area_id)
    end

    # start iteration
    while iteration <= max_iteration && !global_flag_convergence

        Distributed.@everywhere keys(worker_area) begin
            for area in area_id
                # solve local problem and update solution
                result = solve_pmada_model(data_local[area], model_type, optimizer, dopf_method.build_method, solution_processors=dopf_method.post_processors)
                update_data!(data_local[area], result["solution"])
        
                # send data to neighboring areas
                for neighbor in data_local[area]["neighbors"] 
                    shared_data = prepare_shared_data(data_local[area], neighbor)
                    put!(comms[area][neighbor], shared_data)
                end
            end
        end

        Distributed.@everywhere keys(worker_area) begin
            for area in area_id
                # receive data to neighboring areas
                for neighbor in data_local[area]["neighbors"] 
                    received_data = take!(comms[neighbor][area])
                    receive_shared_data!(data_local[area], received_data, neighbor)
                end

                # calculate and share mismatches
                dopf_method.update_method(data_local[area])
                counters = Dict("option"=> data_local[area]["option"], "counter" => data_local[area]["counter"], "mismatch" => data_local[area]["mismatch"])
                if data_local[area]["option"]["termination_measure"] in ["dual_residual", "mismatch_dual_residual"]
                    counters["dual_residual"] = data_local[area]["dual_residual"]
                end
                put!(comms[area][0], deepcopy(counters))
            end
        end

        # receive the mismatches from areas
        for area in areas_id
            counters = take!(comms[area][0])
            global_counters[area] = counters
        end

        # print progress
        print_iteration(global_counters, print_level)

        # update flag convergence and iteration number
        global_flag_convergence = update_global_flag_convergence(global_counters)
        iteration += 1

    end

    # receive the final solution
    Distributed.@everywhere keys(worker_area) begin
        for area in area_id
            # send the area data
            put!(comms[area][0], data_local[area])
        end
    end

    for area in areas_id
        data_area[area] = take!(comms[area][0])
    end

    # close the communication channels
    for i in keys(comms)
        for j in keys(comms[i])
            close(comms[i][j])
        end
    end

    print_convergence(data_area, print_level)
    return data_area
end

"""
    solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using distributed algorithm with central coordinator.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize\\_method_local::Function : initialize the local algorithm parameters and shared variables
    - initialize\\_method_coordinator::Function : initialize the coordinator algorithm parameters and shared variables
    - update\\_method_local::Function : update the local data after each iteration
    - update\\_method_coordinator::Function : update the coordinator data after each iteration
    - build\\_method_local::Function : local problem formulation
    - build\\_method_coordinator::Function : coordinator problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- multiprocessors::Bool=false : enable multiprocessors using available workers. Multiprocessors feature requires loading the PowerModelsADA and the optimizer packages on all the processors using @everywhere using <package_name>.
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous\\_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    # arrange and get areas id
    arrange_areas_id!(data)
    areas_id = get_areas_id(data)

    if length(areas_id) < 2
        error("Number of areas is less than 2, at least 2 areas is needed")
    end
    
    # decompose the system into subsystems
    data_area = Dict{Int64, Any}(0 => decompose_coordinator(data))
    for area in areas_id
        data_area[area] = decompose_system(data, area)
    end

    solve_dopf_coordinated(data_area, model_type, optimizer, dopf_method; print_level=print_level, multiprocessors=multiprocessors, kwargs...)
end

function solve_dopf_coordinated(data::String, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    data = parse_file(data)
    solve_dopf_coordinated(data, model_type, optimizer, dopf_method; print_level=print_level, multiprocessors=multiprocessors, kwargs...)
end

function solve_dopf_coordinated(data::Dict{Int64, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, kwargs...)
    if !multiprocessors
        solve_dopf_coordinated_sp(data, model_type, optimizer, dopf_method; print_level=print_level, kwargs...)
    else
        solve_dopf_coordinated_mp(data, model_type, optimizer, dopf_method; print_level=print_level, kwargs...)
    end
end

"""
    solve_dopf_coordinated_sp(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using distributed algorithm with central coordinator on single-processors.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize\\_method_local::Function : initialize the local algorithm parameters and shared variables
    - initialize\\_method_coordinator::Function : initialize the coordinator algorithm parameters and shared variables
    - update\\_method_local::Function : update the local data after each iteration
    - update\\_method_coordinator::Function : update the coordinator data after each iteration
    - build\\_method_local::Function : local problem formulation
    - build\\_method_coordinator::Function : coordinator problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous\\_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf_coordinated_sp(data_area::Dict{Int64, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)

    # get areas ids
    areas_id = get_areas_id(data_area)
    deleteat!(areas_id, areas_id .== 0)

    # initilize distributed power model parameters
    dopf_method.initialize_method_coordinator(data_area[0], model_type; kwargs...)
    for area in areas_id
        dopf_method.initialize_method_local(data_area[area], model_type; kwargs...)
    end

    ## get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    # initialize the algorithms global counters
    iteration = 0
    flag_convergence = false

    # start iteration
    while iteration <= max_iteration && !flag_convergence

        # solve local area problems in parallel
        info1 = @capture_out begin
            Threads.@threads for area in areas_id
                result = solve_pmada_model(data_area[area], model_type, optimizer, dopf_method.build_method_local, solution_processors=dopf_method.post_processors_local)
                update_data!(data_area[area], result["solution"])
            end
        end

        # share solution of local areas with the coordinator
        for area in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_area[area], 0, serialize = false)
            receive_shared_data!(data_area[0], deepcopy(shared_data), area)
        end

        # solve coordinator problem 
        info2 = @capture_out begin
            result = solve_pmada_model(data_area[0], model_type, optimizer, dopf_method.build_method_coordinator, solution_processors=dopf_method.post_processors_coordinator)
            update_data!(data_area[0], result["solution"])
        end

        # share coordinator solution with local areas
        for area in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_area[0], area, serialize = false)
            receive_shared_data!(data_area[area], deepcopy(shared_data), 0)
        end

        # update local areas and coordinator problems after
        dopf_method.update_method_coordinator(data_area[0])
        for area in areas_id
            dopf_method.update_method_local(data_area[area])
        end

        # print solution
        print_iteration(data_area, print_level, [info1; info2])

        # check global convergence and update iteration counters
        flag_convergence = data_area[0]["counter"]["flag_convergence"]
        iteration += 1

    end

    print_convergence(data_area, print_level)

    return data_area
end

"""
    solve_dopf_coordinated_mp(data::Dict{String, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, multiprocessors::Bool=false, mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, save_data::Vector{String}=[], kwargs...)

Solve OPF problem using distributed algorithm with central coordinator on multiprocessors. Multiprocessors feature requires loading the PowerModelsADA and the optimizer packages on all the processors using @everywhere using <package_name>.
# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- dopf_method::Module : module contains the distributed algorithm methods as follows:
    - initialize\\_method_local::Function : initialize the local algorithm parameters and shared variables
    - initialize\\_method_coordinator::Function : initialize the coordinator algorithm parameters and shared variables
    - update\\_method_local::Function : update the local data after each iteration
    - update\\_method_coordinator::Function : update the coordinator data after each iteration
    - build\\_method_local::Function : local problem formulation
    - build\\_method_coordinator::Function : coordinator problem formulation
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- save_data::Vector{String}=[] : vector contains the keys of the dictionaries to be saved at each iteration in "previous\\_solution". For example, save_data=["solution", "shared_variable", "mismatch"]
- kwargs = includes algorithm-specific and initialization parameters
"""
function solve_dopf_coordinated_mp(data_area::Dict{Int, <:Any}, model_type::DataType, optimizer, dopf_method::Module; print_level::Int64=1, kwargs...)

    # lookup dictionaries for worker-area pairs
    areas_id = get_areas_id(data_area)
    deleteat!(areas_id, areas_id .== 0)
    worker_id = Distributed.workers()
    number_workers = length(worker_id)

    k = 1
    area_worker = Dict()
    for i in areas_id
        if k > number_workers
            k = 1
        end
        area_worker[i] = worker_id[k]
        k += 1 
    end
    worker_area = Dict([i => findall(x -> x==i, area_worker) for i in worker_id if i in values(area_worker)])

    # initiaiate communication channels 
    comms = Dict(0 => Dict(area => Distributed.RemoteChannel(1) for area in areas_id))
    for area1 in areas_id
        comms[area1] = Dict()
        for area2 in [0; areas_id]
            if area1 != area2 
                comms[area1][area2] = Distributed.RemoteChannel(area_worker[area1])
            end
        end
    end

    # initilize distributed power model parameters
    dopf_method.initialize_method_coordinator(data_area[0], model_type; kwargs...)
    for area in areas_id
        dopf_method.initialize_method_local(data_area[area], model_type; kwargs...)
        put!(comms[0][area], data_area[area])
    end

    ## get global parameters
    max_iteration = get(kwargs, :max_iteration, 1000)

    # initialize the algorithms global counters
    iteration = 0
    flag_convergence = false

    # share global variables
    Distributed.@everywhere keys(worker_area) begin
        comms = $comms
        areas_id = $areas_id
        worker_area = $worker_area
        area_worker = $area_worker
        dopf_method = $dopf_method
        model_type = $model_type
        optimizer = $optimizer
        area_id = worker_area[myid()]
        data_local = Dict(area => take!(comms[0][area]) for area in area_id)
    end

    # start iteration
    while iteration <= max_iteration && !flag_convergence
        Distributed.@everywhere keys(worker_area) begin
            for area in area_id
                # solve local problem and update solution
                result = solve_pmada_model(data_local[area], model_type, optimizer, dopf_method.build_method_local, solution_processors=dopf_method.post_processors_local)
                update_data!(data_local[area], result["solution"])
                # send data to coordinator
                shared_data = prepare_shared_data(data_local[area], 0)
                put!(comms[area][0], shared_data)
            end
        end

        # share solution of local areas with the coordinator
        for area in areas_id
            received_data = take!(comms[area][0])
            receive_shared_data!(data_area[0], received_data, area)
        end

        # solve coordinator problem 
        result = solve_pmada_model(data_area[0], model_type, optimizer, dopf_method.build_method_coordinator, solution_processors=dopf_method.post_processors_coordinator)
        update_data!(data_area[0], result["solution"])
        dopf_method.update_method_coordinator(data_area[0])
        # share coordinator solution with local areas
        for area in areas_id
            shared_data = prepare_shared_data(data_area[0], area)
            put!(comms[0][area], shared_data)
        end

        Distributed.@everywhere keys(worker_area) begin
            for area in area_id
                # receive data to neighboring areas
                received_data = take!(comms[0][area])
                receive_shared_data!(data_local[area], received_data, 0)

                # calculate mismatches and update convergence flags
                dopf_method.update_method_local(data_local[area])
            end
        end

        # print solution
        print_iteration_coordinator(data_area, print_level, [])

        # check global convergence and update iteration counters
        flag_convergence = data_area[0]["counter"]["flag_convergence"]
        iteration += 1

    end

    if number_workers > 1
        Distributed.@everywhere keys(worker_area) begin
            for area in area_id
                # send the area data
                put!(comms[area][0], data_local[area])
            end
        end

        for area in areas_id
            data_area[area] = take!(comms[area][0])
        end
    end

    for i in keys(comms)
        for j in keys(comms[i])
            close(comms[i][j])
        end
    end

    print_convergence(data_area, print_level)
    return data_area
end

"initialize dopf parameters"
function initialize_dopf!(data::Dict{String, <:Any}, model_type::DataType; kwargs...)
    # options
    data["option"] = Dict{String, Any}()
    data["option"]["tol"] = get(kwargs, :tol, 1e-4)
    data["option"]["max_iteration"] = get(kwargs, :max_iteration, 1000)
    data["option"]["mismatch_method"] = get(kwargs, :mismatch_method, "norm")
    data["option"]["model_type"] = model_type
    data["option"]["termination_method"] = get(kwargs, :termination_method, "global")
    data["option"]["termination_measure"] = get(kwargs, :termination_measure, "mismatch")

    # counters
    data["counter"] = Dict{String, Any}()
    data["counter"]["iteration"] = Int64(1)
    data["counter"]["flag_convergence"] = false
    data["counter"]["convergence_iteration"] = Int64(0)

    areas_id = get_areas_id(data)
    area_id = get_area_id(data)

    # mismatch 
    data["mismatch"] = Dict{String, Any}()
    data["neighbors"] = [area for area in areas_id if area != area_id]

    # distributed termination method
    if data["option"]["termination_method"] in ["local", "distributed"]
        areas_id = string.(areas_id)
        area_id = string(area_id)
        deleteat!(areas_id, areas_id .== area_id)
        all_areas = string.(get(kwargs, :all_areas, []))
        data["shared_flag_convergence"] = Dict(area => Dict(area_ => false for area_ in all_areas) for area in areas_id)
        data["received_flag_convergence"] = Dict(area => Dict(area_ => false for area_ in all_areas) for area in areas_id)
        data["shared_convergence_iteration"] = Dict(area => 0 for area in areas_id)
        data["received_convergence_iteration"] = Dict(area => 0 for area in areas_id)
        data["counter"]["local_flag_convergence"] = false
        data["option"]["diameter"] = get(kwargs, :diameter, size(all_areas)[1])
    end

    # last solution 
    initialization_method = get(kwargs, :initialization_method, "flat")
    data["solution"] = initialize_all_variable(data, model_type, initialization_method)

    # previous solutions 
    save_data = get(kwargs, :save_data, [])
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
    calc_mismatch!(data::Dict{String, <:Any}; central::Bool=false)
calculate the mismatch and return the area data dictionary with the mismatch as seen by the area. Set central=true if the algorithm uses the optimality condition of a central coordinator.
"""
function calc_mismatch!(data::Dict{String, <:Any}; central::Bool=false)
    area_id = string(get_area_id(data)) 
    mismatch_method = data["option"]["mismatch_method"]
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]

    mismatch = Dict{String, Any}([
        area => Dict{String, Any}([
            variable => Dict{String, Any}([
                idx => central ? (shared_variable_local[area][variable][idx] - (shared_variable_received[area][variable][idx] +shared_variable_local[area][variable][idx] )/2) : (shared_variable_local[area][variable][idx] - shared_variable_received[area][variable][idx])
            for idx in keys(shared_variable_local[area][variable])])
        for variable in keys(shared_variable_local[area])])
    for area in keys(shared_variable_local) if area != area_id && area in keys(shared_variable_received) ])

    if mismatch_method == "norm"
        mismatch[area_id] = LinearAlgebra.norm([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]], 2)
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        mismatch[area_id] = LinearAlgebra.maximum([abs(value) for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]])
    end

    data["mismatch"] = mismatch
end

"""
    calc_dual_residual!(data::Dict{String, <:Any}; central::Bool=false)
calculate the dual redidual as seen by the area. Set central=true if the algorithm uses the optimality condition of a central coordinator.
"""
function calc_dual_residual!(data::Dict{String, <:Any}; central::Bool=false)
    area_id = string(get_area_id(data))
    mismatch_method = data["option"]["mismatch_method"]
    alpha = data["parameter"]["alpha"]
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]

    if data["counter"]["iteration"] == 1
        dual_dual_residual = Dict{String, Any}([
            area => Dict{String, Any}([
                variable => Dict{String, Any}([
                    idx => central ? -alpha* (shared_variable_local[area][variable][idx]+shared_variable_received[area][variable][idx])/2 : -alpha* shared_variable_local[area][variable][idx]
                for idx in keys(shared_variable_local[area][variable])])
            for variable in keys(shared_variable_local[area])])
        for area in keys(shared_variable_local)])
    else
        previous_shared_variable_local = data["previous_solution"]["shared_variable"][end]
        previous_shared_variable_received = data["previous_solution"]["received_variable"][end]
        dual_dual_residual = Dict{String, Any}([
            area => Dict{String, Any}([
                variable => Dict{String, Any}([
                    idx => central ? -alpha * ((shared_variable_local[area][variable][idx]+shared_variable_received[area][variable][idx])/2 - (previous_shared_variable_local[area][variable][idx] +previous_shared_variable_received[area][variable][idx] )/2) : -alpha * (shared_variable_local[area][variable][idx] - previous_shared_variable_local[area][variable][idx])
                for idx in keys(shared_variable_local[area][variable])])
            for variable in keys(shared_variable_local[area])])
        for area in keys(shared_variable_local) ])
    end

    if mismatch_method == "norm"
        dual_dual_residual[area_id] = LinearAlgebra.norm([value for area in keys(dual_dual_residual) if area != area_id for variable in keys(dual_dual_residual[area]) for (idx,value) in dual_dual_residual[area][variable]])
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        dual_dual_residual[area_id] = LinearAlgebra.maximum([abs(value) for area in keys(dual_dual_residual) if area != area_id for variable in keys(dual_dual_residual[area]) for (idx,value) in dual_dual_residual[area][variable]])
    end

    data["dual_residual"] = dual_dual_residual

end

# "get the parameter for each consistency variable as a dictionary"
# function get_parameter(data::Dict{String, <:Any}, parameter::String)
#     if haskey(data, parameter)
#         return data[parameter]
#     else
#         return Dict{String, Any}([area => Dict{String, Any}([variable => Dict{String, Any}([idx => data["parameter"][parameter] for idx in keys(data[parameter][area][variable])]) for variable in keys(data[parameter][area])]) for area in keys(data[parameter])])
#     end
# end

"check flag convergance using mismatch and dual residual"
function flag_convergance(data::Dict{String, <:Any})

    area_id = string(data["area"])
    termination_measure = data["option"]["termination_measure"]

    tol = data["option"]["tol"]
    mismatch = data["mismatch"][area_id]

    if termination_measure in ["dual_residual", "mismatch_dual_residual"]
        tol_dual = data["option"]["tol_dual"]
        dual_residual = data["dual_residual"][area_id]
        flag_convergence = (mismatch < tol && dual_residual < tol_dual)
    else
        flag_convergence = mismatch < tol
    end

    return flag_convergence
end

"check the shared variables of a local area are within tol"
function update_flag_convergence!(data::Dict{String, <:Any})
    area_id = string(data["area"])
    areas_id = string.(get_areas_id(data))
    deleteat!(areas_id, areas_id .== area_id)

    iteration = data["counter"]["iteration"]
    flag_convergence = flag_convergance(data)

    if data["option"]["termination_method"] == "global"
        # the convergence flag is communicated globally
        if flag_convergence && !data["counter"]["flag_convergence"]
            data["counter"]["convergence_iteration"] = iteration
        end
        data["counter"]["flag_convergence"] = flag_convergence
    else # the convergence flag is decided locally
        # Rule 1
        if flag_convergence && !data["counter"]["local_flag_convergence"]
            data["counter"]["convergence_iteration"] = iteration
            data["counter"]["local_flag_convergence"] = flag_convergence

            for area in keys(data["shared_flag_convergence"])
                data["shared_flag_convergence"][area][area_id] = flag_convergence
            end

        end

        # Rule 2
        all_areas = string.(collect(keys(data["shared_flag_convergence"][areas_id[1]])))
        shared_convergence_iteration = maximum([data["counter"]["convergence_iteration"] ; [data["received_convergence_iteration"][area] for area in areas_id] ])

        for area1 in areas_id
            data["shared_convergence_iteration"][area1] = shared_convergence_iteration
            for area2 in all_areas
                if area2 != area_id
                    data["shared_flag_convergence"][area1][area2] = reduce( | , [data["received_flag_convergence"][area][area2] for area in areas_id])
                end
            end
        end

        # Rule 3
        global_flag_convergence = reduce( & , [ val for (area, val) in first(data["shared_flag_convergence"])[2]])

        if global_flag_convergence && (shared_convergence_iteration + data["option"]["diameter"] <= iteration)
            data["counter"]["flag_convergence"] = global_flag_convergence
        end
    end

end

# "calculate the global mismatch based on local mismatch"
# function calc_global_mismatch(data_area::Dict{Int, <:Any}) 
#     mismatch_method = first(data_area)[2]["option"]["mismatch_method"]
#     termination_measure = first(data_area)[2]["option"]["termination_measure"]
#     termination_method = first(data_area)[2]["option"]["termination_method"]

#     if termination_method == "global"
#         if mismatch_method == "norm"
#             if termination_measure in ["dual_residual", "mismatch_dual_residual"]
#                 mismatch = LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
#                 dual_residual = LinearAlgebra.norm([data_area[i]["dual_residual"][string(i)] for i in keys(data_area) if i != 0])
#                 return LinearAlgebra.maximum([mismatch, dual_residual])
#             else
#                 return LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
#             end
#         elseif mismatch_method == "max" || mismatch_method == "maximum"
#             if termination_measure in ["dual_residual", "mismatch_dual_residual"]
#                 mismatch = LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
#                 dual_residual = LinearAlgebra.maximum([data_area[i]["dual_residual"][string(i)] for i in keys(data_area) if i != 0])
#                 return LinearAlgebra.maximum([mismatch, dual_residual])
#             else
#                 return LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
#             end
#         end
#     else
#         if termination_measure in ["dual_residual", "mismatch_dual_residual"]
#             return LinearAlgebra.maximum([[data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0]; [data_area[i]["dual_residual"][string(i)] for i in keys(data_area) if i != 0]])
#         else
#             return LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
#         end
#     end
# end

"calculate the global mismatch based on local mismatch"
function calc_global_mismatch(data_area::Dict{Int, <:Any}) 
    mismatch_method = first(data_area)[2]["option"]["mismatch_method"]

    if mismatch_method == "norm"
        return LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        return LinearAlgebra.maximum([data_area[i]["mismatch"][string(i)] for i in keys(data_area) if i != 0])
    end
end

"calculate the global mismatch based on local mismatch"
function calc_global_dual_residual(data_area::Dict{Int, <:Any}) 
    mismatch_method = first(data_area)[2]["option"]["mismatch_method"]

    if mismatch_method == "norm"
        return LinearAlgebra.norm([data_area[i]["dual_residual"][string(i)] for i in keys(data_area) if i != 0])
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        return LinearAlgebra.maximum([data_area[i]["dual_residual"][string(i)] for i in keys(data_area) if i != 0])
    end
end


"check the flag convergence for all areas and return a global variables"
function update_global_flag_convergence(data_area::Dict{Int64, <:Any})
    if first(data_area)[2]["option"]["termination_method"] == "global"
        termination_measure = first(data_area)[2]["option"]["termination_measure"]
        mismatch = calc_global_mismatch(data_area)
        tol = first(data_area)[2]["option"]["tol"]
        if termination_measure in ["dual_residual", "mismatch_dual_residual"]
            tol_dual = first(data_area)[2]["option"]["tol_dual"]
            dual_residual = calc_global_dual_residual(data_area) 
            return mismatch < tol && dual_residual < tol_dual
        else
            mismatch = calc_global_mismatch(data_area)
            tol = first(data_area)[2]["option"]["tol"]
            return mismatch < tol
        end
    else
        return reduce( &, [data_area[i]["counter"]["flag_convergence"] for i in keys(data_area)])
    end
end 

"print iteration information"
function print_iteration(data::Dict{Int64, <:Any}, print_level::Int64, info_list::Vector=[])
    if print_level > 0
        iteration = first(data)[2]["counter"]["iteration"]-1
        mismatch = calc_global_mismatch(data)
        println("Iteration = $iteration, mismatch = $mismatch")

        if print_level > 1
            for info in info_list
                println(info)
            end
        end
    end
end

"print iteration information"
function print_iteration_coordinator(data::Dict{Int64, <:Any}, print_level::Int64, info_list::Vector=[])
    if print_level > 0
        iteration = data[0]["counter"]["iteration"]-1
        mismatch = data[0]["mismatch"]["0"]
        println("Iteration = $iteration, mismatch = $mismatch")

        if print_level > 1
            for info in info_list
                println(info)
            end
        end
    end
end


# function print_iteration(mismatch::Dict, iteration::Int64, print_level::Int64, info_list::Vector=[])
#     if print_level > 0
#         mismatch = LinearAlgebra.norm([mismatch[i] for i in keys(mismatch)])
#         println("Iteration = $iteration, mismatch = $mismatch")

#         if print_level > 1
#             for info in info_list
#                 println(info)
#             end
#         end
#     end
# end

"print final solution status"
function print_convergence(data::Dict, print_level::Int64)
    if print_level > 0
        iteration = first(data)[2]["counter"]["iteration"]-1
        mismatch = calc_global_mismatch(data)
        tol =first(data)[2]["option"]["tol"]
        flag_convergence = update_global_flag_convergence(data)
        if flag_convergence
            println("*******************************************************")
            println("")
            println("Consistency achieved within $tol mismatch tolerance")
            println("Number of iterations = $iteration")
            
            objective = calc_dist_gen_cost(data)
            println("Objective function value = $objective")
            println("")
            println("*******************************************************")
        else
            println("*******************************************************")
            println("")
            println("Consistency did not achieved within $tol mismatch tolerance and $iteration iteration")
            println("Shared variables mismatch = $mismatch")
            println("")
            println("*******************************************************")
        end
    end
end