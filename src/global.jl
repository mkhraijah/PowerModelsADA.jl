## helper methods to wrap the distributed algorithm to run in a global scope

function decompose_system(data::Dict{String, Any})
    areas_id = get_areas_id(data)
    subsystem_data = Dict{Int64, Any}()
    for i in areas_id
        subsystem_data[i] = decompose_system(data, i)
    end
    return subsystem_data
end

##
function instantiate_dpm_model(subsystem_data::Dict{Int64, <:Any}, pf_model ; setting::Dict{String, <:Any}=[])
    subsystem = Dict{Int64, Any}()
    for i in keys(subsystem_data)
        subsystem[i] = instantiate_dpm_model(subsystem_data[i],pf_model, setting = setting)
    end
    return subsystem
end

##
function initialize_dpm!(subsystem::Dict{Int, <:Any}, optimizer, parameters)
    for i in keys(subsystem)
        initialize_dpm!(subsystem[i], optimizer, parameters)
    end
end


##
function share_solution!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem) # sender subsystem
        for j in keys(subsystem) # receiver subsystem
            if i != j && i in keys(subsystem[j].ext[:primal_shared_variable])
                shared_data = send_shared_data(i, j, subsystem[i], serialize = false)
                    ### Communication ####
                receive_shared_data!(i, shared_data, subsystem[j])
            end
        end
    end
end

##
function calc_mismatch!(subsystem::Dict{Int, <:Any}, p::Int64=2 )
    for i in keys(subsystem)
        calc_mismatch!(subsystem[i])
    end
end

##
function update_flag_convergance!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem)
        update_flag_convergance!(subsystem[i])
    end
end

##
function update_iteration!(pms::Dict)
    for i in keys(pms)
        pms[i].ext[:iteration] += 1
    end
end

##
function update_dual_variable!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem)
        update_dual_variable!(subsystem[i])
    end
end

##
function update_objective!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem)
        update_objective!(subsystem[i])
    end
end

##
function solve_local_model!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem)
        solve_local_model!(subsystem[i])
    end
end

##
function update_primal_variable!(subsystem::Dict{Int, <:Any})
    for i in keys(subsystem)
        update_primal_variable!(subsystem[i])
    end
end

##
check_mismatch(subsystem::Dict,p::Int64=2) = sum(subsystem[i].ext[:mismatch][end][i]^p for i in keys(subsystem))^1/p

##
function run_dopf(data::Dict, pf_model, optimizer, setting::Dict{String, <:Any}=Dict{String, <:Any}(), parameters::Real=100; verbose = true)

    ## Check the power flow model
    pf_model= pf_formulation(pf_model)
    ## Decompose the system into several subsystem return objects of power model subtype
    subsystem_data = decompose_system(data)
    subsystem = instantiate_dpm_model(subsystem_data , pf_model, setting = setting)
    initialize_dpm!(subsystem, optimizer, parameters)

    ## Initialaize the algorithms counters
    iteration = 1
    flag_convergance = false
    max_iteration = get(setting,"max_iteration",1000)
    areas_id = get_areas_id(data)

    ## start iteration
    while iteration < max_iteration && flag_convergance == false

        update_dual_variable!(subsystem)
        update_objective!(subsystem)
        solve_local_model!(subsystem)
        # Check if there is a solution and terminiate if there is no solution
        status = check_solution(subsystem)
        have_solution = reduce( & , [status["solution_status"][i] for i in keys(subsystem)] )
        if !have_solution
            failed_subsystems = [ i for (i,st) in status["solution_status"] if !st]
            println("Subsystem/s $failed_subsystems failed to obtain a solution due to $([status["subsystem_termination"][i] for i in failed_subsystems]).")
            break
        end

        update_primal_variable!(subsystem)
        share_solution!(subsystem)
        calc_mismatch!(subsystem,2)
        update_flag_convergance!(subsystem)

        flag_convergance = check_flag_convergance(subsystem)

        if verbose
            mismatch = check_mismatch(subsystem)
            println("Iteration = $iteration, mismatch = $mismatch")
            if flag_convergance
                println("Consistency achived within $(setting["tol"])")
            end
        end

        iteration += 1
        update_iteration!(subsystem)

    end

    return subsystem
end

##  function to check the status of the local solutions
function check_solution(subsystem::Dict)
    status = Dict("solution_status" => Dict(), "subsystem_termination" => Dict())
    for (i,system) in subsystem
        status["solution_status"][i] = JuMP.has_values(subsystem[i].model)
        status["subsystem_termination"][i] = JuMP.termination_status(subsystem[i].model)
    end
    return status
end

## Compare the distributed algorithm solutoin with PowerModels centralized solution
function compare_solution(data::Dict, subsystem::Dict, pf_model, optimizer)

    pf_model = pf_formulation(pf_model)
    _PM.standardize_cost_terms!(data, order=2)
    # Get distributed OPF solutoin from agents
    pg_dist = Dict{Any, Any}(j => subsystem[i].solution["gen"][j]["pg"] for i in keys(subsystem) for j in keys(subsystem[i].solution["gen"]))

    # Solve Centralized OPF
    Central_solution = _PM.run_opf(data, pf_model, optimizer)

    # Calculate objective function
    Obj_cent = Central_solution["objective"]
    Obj_dist = sum(gen["cost"][1]*pg_dist[i]^2 + gen["cost"][2]*pg_dist[i]+ gen["cost"][3] for (i,gen) in data["gen"] if gen["gen_status"] == 1)


    # Calculate optimility gap
    Optimility_Gap = (Obj_dist - Obj_cent)/ Obj_cent * 100
    println("")
    println("Relative Error with PM = $Optimility_Gap")
    println("")
    return Optimility_Gap
end
