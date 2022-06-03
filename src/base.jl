###############################################################################
#              Base method for all distirbuted OPF algorithms                 #
###############################################################################

## method to initialize dpm parameters and optimizer
function initialize_dpm!(data::Dict{String, <:Any}, model_type::Type; alpha::Real, tol::Float64=1e-4, max_iteration::Int64=1000)

    # initiate primal and dual shared variables
    variable_shared(data, model_type)

    # initiate distributed algorithm parameters
    data["iteration"] = Int64(1)
    data["flag_convergance"] = false
    data["mismatch"] = Vector{Dict{String, Any}}()
    data["alpha"] = alpha
    data["tol"] = tol
    data["max_iteration"] = max_iteration

end


# ## wraping PowerModel.instantiate_model with dopf setting
# function update_subproblem(data::Dict{String, <:Any}, model_type, build_method::Function; alpha::Real=1000, tol::Float64=1e-4, max_iteration::Int64=1000)
#
#     ## check the power flow model
#     model_type= pf_formulation(model_type)
#
#     # ensure the subsystem data includes the area id
#     if !haskey(data,"area")
#         error("No area id is provided in the data")
#     end
#
#     setting = Dict{String,Any}("alpha" => alpha, "tol" => tol, "max_iteration" => max_iteration)
#
#
#     _PM.instantiate_model(data, model_type, build_method, setting = setting)
# end


# ## solve local dopf called every iteration
# function solve_subproblem!(pm::AbstractPowerModel, optimizer)
#     JuMP.set_optimizer(pm.model, optimizer)
#     JuMP.set_silent(pm.model)
#     optimize_model!(pm)
#     update_primal_variable!(pm)
#     update_data!(pm.data, pm.solution)
# end


# ## wrapping for updata_data! in PowerModels
# function update_solution!(data::Dict{String,<:Any}, pm::AbstractPowerModel)
#     update_data!(data, pm.solution)
#     data["shared_primal"] = pm.data["shared_primal"]
#     data["shared_dual"] = pm.data["shared_dual"]
#     update_iteration!(data)
# end

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

    mismatch[area_id] = norm([value for area in keys(mismatch) if area != area_id for comp in keys(mismatch[area]) for ids in keys(mismatch[area][comp]) for (vstring,value) in mismatch[area][comp][ids]],p)

    push!(data["mismatch"], mismatch)
end

## Check the shared variables of a local area are within tol
function update_flag_convergance!(data::Dict{String, <:Any}, tol::Float64)
    area_id = string(data["area"])
    mismatch = data["mismatch"][end][area_id]
    data["flag_convergance"] = mismatch < tol
end
