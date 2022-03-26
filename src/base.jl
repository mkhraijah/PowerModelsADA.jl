

## wrapper for PowerModel.instantiate_model with dopf setting
function instantiate_dpm_model(data::Dict{String, Any}, pf_model; setting::Dict=Dict())

    if !haskey(data,"area")
        error("No area id is provided in the data")
    end
    if isempty(setting)
        setting = set_setting()
    end
    pf_model = DPM.pf_formulation(pf_model)

    _PM.instantiate_model(data, pf_model, build_dopf, setting = setting)
end

function instantiate_dpm_model(data::Dict{String, Any}, area::Int64, pf_model; setting::Dict{String, Any}=Dict{String, Any}())
    data_area = decompose_system(data,area)
    instantiate_dpm_model(data_area, setting=setting)
end

## build method for Distributed PowerModel
function build_dopf(pm::AbstractPowerModel)

    # define variables
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

    # initiate primal and dual shared variables in pm.ext
    variable_shared(pm)

    # define constraints
    _PM.constraint_model_voltage(pm)

    for i in ids(pm, :ref_buses)
        _PM.constraint_theta_ref(pm, i)
    end

    for i in ids(pm, :bus)
        _PM.constraint_power_balance(pm, i)
    end

    for i in ids(pm, :branch)
        _PM.constraint_ohms_yt_from(pm, i)
        _PM.constraint_ohms_yt_to(pm, i)
        _PM.constraint_thermal_limit_from(pm, i)
        _PM.constraint_thermal_limit_to(pm, i)
        _PM.constraint_voltage_angle_difference(pm, i)
    end

    for i in ids(pm, :dcline)
        _PM.constraint_dcline_power_losses(pm, i)
    end

end

""" I'm not sure if this is the best method to wrap the setting into a dict """
## method to wrap dpm settings into a dictionary
function set_setting(distributed_algorithm::String="APP", tol::Float64= 1e-4, max_iteration::Int64=1000)

    setting = Dict{String,Any}("distributed_algorithm" => distributed_algorithm, "tol" => tol, "max_iteration" => max_iteration)

end

## method to initialize dpm parameters and optimizer
function initialize_dpm!(pm::AbstractPowerModel, optimizer , parameters...)
    pm.ext[:iteration] = Int64(0)
    pm.ext[:flag_convergance] = false
    pm.ext[:mismatch] = Vector{Dict{Int64, Any}}()
    pm.ext[:parameter] = parameters
    pm.ext[:optimizer] = optimizer

    JuMP.set_optimizer(pm.model, optimizer)
    JuMP.set_silent(pm.model)

    if pm.setting["distributed_algorithm"] == "ATC"
        pm.ext[:beta] = 1
    end
end


## solve local dopf called every iteration
function solve_local_model!(pm::AbstractPowerModel)
    _IM.optimize_model!(pm)
end

## update iteration
function update_iteration!(pm::AbstractPowerModel)
    pm.ext[:iteration] += 1
end

## calculate the mismatch and store it in pm.ext[:mismatch]
# for iteration k area i, pm.ext[:mismatch][k][j] is
# if i = j, the p-norm of all mismatches
# if i ≆ j, the actual mismatches for every shared variable
function calc_mismatch!(pm::AbstractPowerModel, p::Int64=2 )
    area = pm.data["area"]
    primal_variable = pm.ext[:primal_shared_variable]
    mismatch = Dict{Int64, Any}([i => Dict{Symbol, Any}([j =>Dict{Any, Any}([ k => primal_variable[area][j][k] - primal_variable[i][j][k] for k in keys(primal_variable[i][j])]) for j in keys(primal_variable[i])]) for i in keys(primal_variable) if i != area ])
    mismatch[area] = norm([value for i in keys(mismatch) for j in keys(mismatch[i]) for (k,value) in mismatch[i][j]],p)
    append!(pm.ext[:mismatch], [mismatch])
end

## Check the shared variables of a local area are within tol
function update_flag_convergance!(pm::AbstractPowerModel)
    area = pm.data["area"]
    tol = pm.setting["tol"]
    mismatch = pm.ext[:mismatch][end][area]
    pm.ext[:flag_convergance] = mismatch < tol
end

##
function check_flag_convergance(pms::Dict{Int64, Any})
    flag_convergance = reduce( & , [pms[i].ext[:flag_convergance] for i in keys(pms)])
    return flag_convergance
end
