###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

## method to initialize dpm parameters and optimizer
function initialize_dpm!(data::Dict{String, <:Any}, pf_model)

    # initiate primal and dual shared variables
    variable_shared(data, pf_model)

    # initiate distributed algorithm parameters
    data["iteration"] = Int64(1)
    data["flag_convergance"] = false
    data["mismatch"] = Vector{Dict{Int64, Any}}()

end


## wraping PowerModel.instantiate_model with dopf setting
function update_subproblem(data::Dict{String, <:Any}, pf_model, build_method::Function; alpha::Real=1000, tol::Float64=1e-4, max_iteration::Int64=1000)

    ## check the power flow model
    pf_model= pf_formulation(pf_model)

    # ensure the subsystem data includes the area id
    if !haskey(data,"area")
        error("No area id is provided in the data")
    end

    setting = Dict{String,Any}("alpha" => alpha, "tol" => tol, "max_iteration" => max_iteration)


    _PM.instantiate_model(data, pf_model, build_method, setting = setting)
end


## solve local dopf called every iteration
function solve_subproblem!(pm::AbstractPowerModel, optimizer)
    JuMP.set_optimizer(pm.model, optimizer)
    JuMP.set_silent(pm.model)
    optimize_model!(pm)
    update_primal_variable!(pm)
    update_data!(pm.data, pm.solution)
end


## wrapping for updata_data! in PowerModels
function update_solution!(data::Dict{String,<:Any}, pm::AbstractPowerModel)
    update_data!(data, pm.solution)
    data["shared_primal"] = pm.data["shared_primal"]
    data["shared_dual"] = pm.data["shared_dual"]
    update_iteration!(data)
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
    area = data["area"]
    primal_variable = data["shared_primal"]
    mismatch = Dict{Int64, Any}([i => Dict{Symbol, Any}([j =>Dict{Any, Any}([ k => primal_variable[area][j][k] - primal_variable[i][j][k] for k in keys(primal_variable[i][j])]) for j in keys(primal_variable[i])]) for i in keys(primal_variable) if i != area ])
    mismatch[area] = norm([value for i in keys(mismatch) for j in keys(mismatch[i]) for (k,value) in mismatch[i][j]],p)
    push!(data["mismatch"], mismatch)
end

## Check the shared variables of a local area are within tol
function update_flag_convergance!(data::Dict{String, <:Any}, tol::Float64)
    area = data["area"]
    mismatch = data["mismatch"][end][area]
    data["flag_convergance"] = mismatch < tol
end
