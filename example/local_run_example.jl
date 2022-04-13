
path = joinpath(@__DIR__,"..")
cd(path)
using Pkg
Pkg.activate(".")
using Ipopt

## Call Packages
include("../src/DistributedPowerModels.jl")
DPM = DistributedPowerModels
## Read case with partition file and return dictionary of the paritioned case
case_path = "test/data/case300.m"
data = DPM.parse_file(case_path)

partition_path= "test/data/case300_3areas.csv"
DPM.assign_area!(data, partition_path)

## Settings
alpha = 10000
max_iteration = 10000
tol = 1e-4
distributed_algorithm = "APP" # Select algorithm among APP ADMM ATC
pf_model = "DC" # Select model among DC AC SOC SDP QC

settings = DPM.set_setting(distributed_algorithm, tol, max_iteration)

## Decompose the system into several subsystem return PowerModel

areas_id = DPM.get_areas_id(data) # return a vector with areas ids
subsystem_data = Dict{Int,Any}()
subsystem = Dict{Int,Any}()
optimizer = Dict{Int,Any}()

for i in areas_id
    optimizer[i] = Ipopt.Optimizer
    subsystem_data[i] = DPM.decompose_system(data, i)
    subsystem[i] = DPM.instantiate_dpm_model(subsystem_data[i], pf_model, setting = settings)
    DPM.initialize_dpm!(subsystem[i], optimizer[i], alpha)
end

## Initialaize the algorithms counters
iteration = 1
flag_convergance = false

## Start the algorithm iteration

@time while iteration < max_iteration && flag_convergance == false
    ## Solve local problem and update shared variables
    for i in areas_id
        DPM.update_dual_variable!(subsystem[i])
        DPM.update_objective!(subsystem[i])
        DPM.solve_local_model!(subsystem[i])
        DPM.update_primal_variable!(subsystem[i])
    end

    ## Share solution
    for i in areas_id # sender subsystem
        for j in areas_id # receiver subsystem
            if i != j && i in keys(subsystem[j].ext[:primal_shared_variable])
                shared_data = DPM.send_shared_data(i, j, subsystem[i], serialize = false)
                    ### Communication ####
                DPM.receive_shared_data!(i, shared_data, subsystem[j])
            end
        end
    end

    ## Calculate mismatches and update convergance flags
    for i in areas_id
        DPM.calc_mismatch!(subsystem[i],2)
        DPM.update_flag_convergance!(subsystem[i])
    end

    ## Calculate global mismatch
    mismatch =  DPM.check_mismatch(subsystem)
    println("Iteration Number = $iteration, mismatch = $mismatch")

    ## Check global convergance and update iteration counters
    flag_convergance = DPM.check_flag_convergance(subsystem)
    iteration += 1
    DPM.update_iteration!(subsystem)

end

DPM.compare_solution(data, subsystem, pf_model, optimizer[1])
