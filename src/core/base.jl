###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

"""
    solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, initialize_method::Function, update_method::Function, build_method::Function ; communication_noise::Function=add_communication_noise!, tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, kwargs...)

Solve the distributed OPF problem. The distributed algorithm is defined by the build_method and update_method.
"""

function solve_dopf(data::Dict{String, <:Any}, model_type::DataType, optimizer, initialize_method::Function, update_method::Function, build_method::Function ; communication_noise::Function=add_communication_noise!, tol::Float64=1e-4, max_iteration::Int64=1000, verbose::Bool=true, kwargs...)

    ## obtain areas idx
    areas_id = get_areas_id(data)

    ## decompose the system into subsystems
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    ## initilize distributed power model parameters
    for i in areas_id
        initialize_method(data_area[i], model_type; tol=tol, max_iteration=max_iteration, kwargs...)
    end

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        ## solve local problem and update solution
        for i in areas_id
            update_method(data_area[i])
            result = _PM.solve_model(data_area[i], model_type, optimizer, build_method)
            update_solution!(data_area[i], result["solution"], model_type)
        end

        ## share solution and add noise to the shared data
        for i in areas_id # sender subsystem
            for j in areas_id # receiver subsystem
                if i != j && string(i) in keys(data_area[j]["shared_primal"])
                    shared_data = send_shared_data(i, j, data_area[i], serialize = false)
                    communication_noise(shared_data)
                    receive_shared_data!(i, shared_data, data_area[j])
                end
            end
        end

        ## calculate mismatches and update convergance flags
        for i in areas_id
            calc_mismatch!(data_area[i])
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
    decompose_system(data::Dict{String, <:Any})

Method to decompose a system into subsystem defined by bus area.
"""
function decompose_system(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end
    return data_area
end

"initialize dopf parameters"
function initialize_dopf!(data::Dict{String, <:Any}, model_type::DataType;  tol::Float64=1e-4, max_iteration::Int64=1000)
    # initiate primal and dual shared variables
    variable_shared(data, model_type)

    # initiate distributed algorithm parameters
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

"communication noise model"
function add_communication_noise!(shared_data) 
    # a placeholder to add communicaiton noise modles 
end


"""
    calc_mismatch!(data::Dict{String, <:Any}, p::Int64=2)
calculate the mismatch using p-norm and return the area data dictionary with the mismatch as seen by the area. The mismatches are stored in a vector indexed by the iteraiton number.
"""
function calc_mismatch!(data::Dict{String, <:Any}, method::String="norm"; p::Int64=2 )
    area_id = string(data["area"])
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
    elseif method == "max" || "maximum"
        mismatch[area_id] = LinearAlgebra.maximum([value for area in keys(mismatch) if area != area_id for variable in keys(mismatch[area]) for (idx,value) in mismatch[area][variable]])
    end

    data["mismatch"] = mismatch
end

"check the shared variables of a local area are within tol"
function update_flag_convergance!(data::Dict{String, <:Any}, tol::Float64)
    area_id = string(data["area"])
    mismatch = data["mismatch"][area_id]
    data["flag_convergance"] = mismatch < tol
end

"update iteration"
function update_iteration!(data::Dict{String, <:Any})
    data["iteration"] += 1
end

"check the flag convergence for all subsystem and return a global variables"
function update_global_flag_convergance(data_area::Dict{Int, <:Any})
    flag_convergance = reduce( & , [data_area[i]["flag_convergance"] for i in keys(data_area)])
    return flag_convergance
end

"calculate the global mismatch based on local mismatch"
calc_global_mismatch(data_area::Dict{Int, <:Any}, p::Int64=2) = LinearAlgebra.norm([data_area[i]["mismatch"][string(i)] for i in keys(data_area)], p)