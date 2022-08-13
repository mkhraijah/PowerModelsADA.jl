###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

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

## method to initialize dopf parameters and optimizer
function initialize_dopf!(data::Dict{String, <:Any}, model_type::Type;  tol::Float64=1e-4, max_iteration::Int64=1000, kwargs...)

    # initiate primal and dual shared variables
    variable_shared(data, model_type)

    # initiate distributed algorithm parameters
    data["iteration"] = Int64(1)
    data["flag_convergance"] = false
    data["mismatch"] = Vector{Dict{String, Any}}()
    data["tol"] = tol
    data["max_iteration"] = max_iteration

end

## update iteration
function update_iteration!(data::Dict{String, <:Any})

    data["iteration"] += 1

end

## calculate the mismatch and store it in data
# for iteration k area i, data["mismatch"][k][j] is
# the p-norm of all mismatches, if i = j
# the actual mismatches for every shared variable otherwise
function calc_mismatch!(data::Dict{String, <:Any}, p::Int64=2 )

    area_id = string(data["area"])
    primal_variable = data["shared_primal"]

    mismatch = Dict{String, Any}([
    area => Dict{String, Any}([
    comp => Dict{String, Any}([
    ids => Dict{String, Any}([
    vstring => primal_variable[area_id][comp][ids][vstring] - primal_variable[area][comp][ids][vstring] for vstring in keys(primal_variable[area][comp][ids])]) for ids in keys(primal_variable[area][comp])]) for comp in keys(primal_variable[area])]) for area in keys(primal_variable) if area != area_id ])

    mismatch[area_id] = LinearAlgebra.norm([value for area in keys(mismatch) if area != area_id for comp in keys(mismatch[area]) for ids in keys(mismatch[area][comp]) for (vstring,value) in mismatch[area][comp][ids]],p)

    push!(data["mismatch"], mismatch)
end

## Check the shared variables of a local area are within tol
function update_flag_convergance!(data::Dict{String, <:Any}, tol::Float64)
    area_id = string(data["area"])
    mismatch = data["mismatch"][end][area_id]
    data["flag_convergance"] = mismatch < tol
end

## Check the flag convergence for all subsystem and return a global variables
function check_flag_convergance(data_area::Dict{Int, <:Any})
    flag_convergance = reduce( & , [data_area[i]["flag_convergance"] for i in keys(data_area)])
    return flag_convergance
end

## Calculate the global mismatch based on local mismatch
calc_global_mismatch(data_area::Dict{Int, <:Any}, p::Int64=2) = LinearAlgebra.norm([data_area[i]["mismatch"][end][string(i)] for i in keys(data_area)], p)


## wrapping method to run distributed algorithms
function run_dopf(data::Dict{String, <:Any}, model_type, build_method::Function, update_method::Function, optimizer; initialize_method::Function=initialize_dopf!, tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, kwargs...)

    ## Obtain areas ids
    areas_id = get_areas_id(data)

    ## Decompose the system into several subsystem return PowerModel
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    ## Initilize distributed power model parameters
    for i in areas_id
        initialize_method(data_area[i], model_type::Type, tol=tol, max_iteration=max_iteration, kwargs=kwargs)
    end

    ## Initialaize the algorithms counters
    iteration = 1
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        ## solve local problem and update solution
        for i in areas_id
            update_method(data_area[i])
            result = _PM.solve_model(data_area[i], model_type, optimizer, build_method)
            _PM.update_data!(data_area[i], result["solution"])
            update_shared_primal!(data_area[i], result["solution"])
        end

        ## Share solution
        for i in areas_id # sender subsystem
            for j in areas_id # receiver subsystem
                if i != j && string(i) in keys(data_area[j]["shared_primal"])
                    shared_data = send_shared_data(i, j, data_area[i], serialize = false)
                        ### Communication ####
                    receive_shared_data!(i, shared_data, data_area[j])
                end
            end
        end

        ## Calculate mismatches and update convergance flags
        for i in areas_id
            calc_mismatch!(data_area[i],2)
            update_flag_convergance!(data_area[i], tol)
        end

        ## Check global convergance and update iteration counters
        flag_convergance = check_flag_convergance(data_area)

        if verbose
            mismatch = calc_global_mismatch(data_area)
            println("Iteration = $iteration, mismatch = $mismatch")
            if flag_convergance
                println("Consistency achived within $tol")
            end
        end

        iteration += 1

    end

    return data_area
end
