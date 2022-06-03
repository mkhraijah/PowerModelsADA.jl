
path = joinpath(@__DIR__,"..")
cd(path)
using Pkg
Pkg.activate(".")
using Ipopt
using JuMP

## Call Packages
include("../src/DistributedPowerModels.jl")
DPM = DistributedPowerModels

## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case14.m"
data = DPM.parse_file(case_path)

# partition_path= "test/data/case300_3areas.csv"
# DPM.assign_area!(data, partition_path)

## Settings
alpha = 1000
max_iteration = 1000
tol = 1e-4
model_type = "DC" # Select model among DC AC SOC SDP QC
optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
build_method = DPM.build_dopf_admm


##  Distributed algorithm settings
areas_id = DPM.get_areas_id(data) # return a vector with areas ids
data_area = Dict{Int,Any}()
model_type = DPM.pf_formulation(model_type)


## Decompose the system into several subsystem return PowerModel
for i in areas_id
    data_area[i] = DPM.decompose_system(data, i)
end

## Initilize distributed power model parameters
for i in areas_id
    DPM.initialize_dpm!(data_area[i], model_type, alpha=1000, tol=1e-4, max_iteration =1000)
end

## Initialaize the algorithms counters
iteration = 1
flag_convergance = false

## Start the algorithm iteration

while iteration < max_iteration && flag_convergance == false

    ## Solve local problem and update shared variables
    for i in areas_id
        pm = DPM.update_subproblem(data_area[i], pf_model, build_method, alpha = alpha, tol = tol, max_iteration = max_iteration)
        DPM.solve_subproblem!(pm, optimizer)
        DPM.update_solution!(data_area[i], pm)
    end

    ## Share solution
    for i in areas_id # sender subsystem
        for j in areas_id # receiver subsystem
            if i != j && i in keys(data_area[j]["shared_primal"])
                shared_data = DPM.send_shared_data(i, j, data_area[i], serialize = false)
                    ### Communication ####
                DPM.receive_shared_data!(i, shared_data, data_area[j])
            end
        end
    end

    ## Calculate mismatches and update convergance flags
    for i in areas_id
        DPM.calc_mismatch!(data_area[i],2)
        DPM.update_flag_convergance!(data_area[i], tol)
    end

    ## Calculate global mismatch
    mismatch =  DPM.calc_global_mismatch(data_area)
    println("Iteration Number = $iteration, mismatch = $mismatch")

    ## Check global convergance and update iteration counters
    flag_convergance = DPM.check_flag_convergance(data_area)
    iteration += 1

end

DPM.compare_solution(data, data_area, pf_model, optimizer)







    ## Obtain areas ids
areas_id = get_areas_id(data)

    ## Decompose the system into several subsystem return PowerModel
data_area = Dict{Int64, Any}()
for i in areas_id
    data_area[i] = decompose_system(data, i)
end

    ## Initilize distributed power model parameters
for i in areas_id
    initialize_method!(data_area[i], model_type, alpha=alpha, tol=tol, max_iteration=max_iteration)
end

    ## Initialaize the algorithms counters
iteration = 1
flag_convergance = false

    ## start iteration
while iteration < max_iteration && flag_convergance == false

        ## solve local problem and update solution
    for i in areas_id
        update_method(data_area[i])
        result = solve_model(data_area[i], model_type, optimizer, build_method, solution_processors=[update_primal_variable!])
        update_data!(data, result["solution"])
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
