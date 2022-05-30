###############################################################################
#     Methods to wrap the distributed algorithm to run in a global scope      #
###############################################################################

function decompose_system(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end
    return data_area
end

##
function initialize_dpm!(data_area::Dict{Int64, <:Any}, optimizer, pf_model)
    for i in keys(data_area)
        initialize_dpm!(data_area[i], optimizer, pf_model)
    end
end

##
function update_subproblem(data_area::Dict{Int64, <:Any}, pf_model, build_method::Function, alpha::Real=1000, tol::Float64=1e-4, max_iteration::Int64=1000)
    pms = Dict{Int64, Any}()
    for i in keys(data_area)
        pms[i] = instantiate_dpm_model(data_area[i], pf_model, build_method, ; alpha = alpha, tol = tol, max_iteration = max_iteration)
    end
    return pms
end




##
function share_solution!(data_area::Dict{Int64, <:Any})
    for i in keys(data_area) # sender subsystem
        for j in keys(data_area) # receiver subsystem
            if i != j && i in keys(data_area[j]["shared_primal"])
                shared_data = send_shared_data(i, j, data_area[i], serialize = false)
                    ### Communication ####
                receive_shared_data!(i, shared_data, data_area[j])
            end
        end
    end
end

##
function calc_mismatch!(data_area::Dict{Int64, <:Any}, p::Int64=2 )
    for i in keys(data_areas)
        calc_mismatch!(data_area[i],p)
    end
end

##
function update_flag_convergance!(data_area::Dict{Int64, <:Any}, tol::Float64)
    for i in keys(data_area)
        update_flag_convergance!(data_area[i], tol)
    end
end

##
function update_iteration!(data_area::Dict{Int64, <:Any})
    for i in keys(data_area)
        data_area[i]["iteration"] += 1
    end
end

##
function solve_local_model!(pms::Dict{Int, <:Any}, optimizer)
    for i in keys(pms)
        solve_subproblem!(pms[i], optimizer)
    end
end


##
calc_global_mismatch(data_area::Dict{Int, <:Any}, p::Int64=2) = norm([data_area[i]["mismatch"][end][i] for i in keys(data_area)])

##
function check_flag_convergance(data_area::Dict{Int, <:Any})
    flag_convergance = reduce( & , [data_area[i]["flag_convergance"] for i in keys(data_area)])
    return flag_convergance
end


## wrapping method to run distributed algorithms
function run_dopf(data::Dict{String, <:Any}, pf_model, build_method::Function, optimizer, alpha::Real=1000; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true)

    ## Check the power flow model
    pf_model= pf_formulation(pf_model)

    ## Obtain areas ids
    areas_id = get_areas_id(data)

    ## Decompose the system into several subsystem return PowerModel
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    ## Initilize distributed power model parameters
    for i in areas_id
        initialize_dpm!(data_area[i], pf_model)
    end

    ## Initialaize the algorithms counters
    iteration = 1
    flag_convergance = false

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        ## solve local problem and update solution
        for i in areas_id
""" simplify things here and make depends on data only (not pm)
make it four steps
step 1 update the dual variable
step 2 optimize power models using the specified build method
step 3 update solution (require converting data into string=>float)
            """
            pm = update_subproblem(data_area[i], pf_model, build_method, alpha = alpha, tol = tol, max_iteration = max_iteration)
            solve_subproblem!(pm, optimizer)
            update_solution!(data_area[i], pm)
        end

        ## Share solution
        for i in keys(data_area) # sender subsystem
            for j in keys(data_area) # receiver subsystem
                if i != j && i in keys(data_area[j]["shared_primal"])
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
""" update the original data file (using update_data! after converting data into string=>float)""" 
    return data_area
end


## Compare the distributed algorithm solutoin with PowerModels centralized solution
function compare_solution(data, data_area, pf_model, optimizer)

    pf_model = pf_formulation(pf_model)

    # Update generator values
    for (i,gen) in data["gen"]
        area = data["bus"]["$(gen["gen_bus"])"]["area"]
        gen["pg"] = data_area[area]["gen"][i]["pg"]
    end

    # Solve Centralized OPF
    Central_solution = _PM.run_opf(data, pf_model, optimizer)

    # Calculate objective function
    Obj_dist = _PM.calc_gen_cost(data)
    Obj_cent = Central_solution["objective"]

    # Calculate optimility gap
    Relative_Error = (Obj_dist - Obj_cent)/ Obj_cent * 100
    return Relative_Error
end
