###############################################################################
#                     Build methods for ATC algorithm                        #
###############################################################################

## solve the distributed OPF problem using ATC algorithm
function run_dopf_atc(data::Dict{String, <:Any}, model_type::Type, optimizer; alpha::Real=1.05, beta::Real=1, tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true)

    run_dopf(data, model_type, build_dopf_atc, update_atc!, optimizer, initialize_method=initialize_dopf_atc! , tol=tol, max_iteration=max_iteration, verbose=verbose, alpha=alpha, beta=beta)

end

## method to inilitlize the ATC algorithm
function initialize_dopf_atc!(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs)

    initialize_dopf!(data, model_type, tol=tol, max_iteration=max_iteration, kwargs=kwargs)
    data["alpha"] = kwargs[:alpha]
    if haskey(kwargs,:beta)
        data["beta"] = kwargs[:beta]
    else
        data["beta"] = 1
    end
end

## build method for Distributed PowerModel using ADMM algorithm
function build_dopf_atc(pm::AbstractPowerModel)

    # define variables
    _PM.variable_bus_voltage(pm)
    _PM.variable_gen_power(pm)
    _PM.variable_branch_power(pm)
    _PM.variable_dcline_power(pm)

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

    objective_min_fuel_and_consensus!(pm, objective_atc!)
end

## method to set the ATC algorithm objective
function objective_atc!(pm::AbstractPowerModel)

    ## ATC parameters
    beta = pm.data["beta"]

    ## data
    area_id = string(get_area_id(pm))
    primal_variable = pm.data["shared_primal"]
    dual_variable = pm.data["shared_dual"]

    ## objective function
    objective = JuMP.objective_function(pm.model)
    for area in keys(primal_variable)
        if area != area_id
            for comp in keys(primal_variable[area])
                for ids in keys(primal_variable[area][comp])
                    for vstring in keys(primal_variable[area][comp][ids])

                        v = pm.sol[:it][:pm][:nw][0][Symbol(comp)][parse(Int64,ids)][Symbol(vstring)]
                        v_central = (primal_variable[area_id][comp][ids][vstring] + primal_variable[area][comp][ids][vstring])/2
                        v_dual = dual_variable[area][comp][ids][vstring]

                        objective += (beta * (v - v_central))^2 + v_dual * (v - v_central)
                    end
                end
            end
        end
    end

    JuMP.@objective(pm.model, Min,  objective)
end


## method to update the dual variable value
function update_atc!(data::Dict{String, <:Any})

    ## ATC parameters
    if !haskey(data,"beta")
        data["beta"] = 1
    end

    alpha = data["alpha"]
    beta  = data["beta"]

    ## data
    area_id = string(get_area_id(data))
    primal_variable = data["shared_primal"]
    dual_variable = data["shared_dual"]

    ## update dual variable
    ## update dual variable
    for area in keys(dual_variable)
        for comp in keys(dual_variable[area])
            for ids in keys(dual_variable[area][comp])
                for vstring in keys(dual_variable[area][comp][ids])

                    v_primal = primal_variable[area_id][comp][ids][vstring]
                    v_central = (primal_variable[area_id][comp][ids][vstring] + primal_variable[area][comp][ids][vstring])/2
                    v_dual = dual_variable[area][comp][ids][vstring]

                    data["shared_dual"][area][comp][ids][vstring] = v_dual  + 2 * beta^2 * (v_primal - v_central)
                end
            end
        end
    end


    ## update ATC parameter
    data["beta"] *= alpha

end
