###############################################################################
#                     Build methods for APP algorithm                        #
###############################################################################

## solve the distributed OPF problem using APP algorithm
function run_dopf_app(data::Dict{String, <:Any}, model_type::Type, optimizer; alpha::Real=1000, beta::Real=0, gamma::Real=0, tol::Float64=1e-4, max_iteration::Int64=1000, verbose=true)

    run_dopf(data, model_type, build_dopf_app, update_app!, optimizer, initialize_method=initialize_dopf_app!, tol=tol, max_iteration=max_iteration, verbose=verbose, alpha=alpha, beta=beta, gamma=gamma)

end

## method to inilitlize the APP algorithm
function initialize_dopf_app!(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs)

    initialize_dopf!(data, model_type, tol=tol, max_iteration=max_iteration, kwargs=kwargs)
    data["alpha"] = kwargs[:alpha]
    data["beta"] = 2*kwargs[:alpha]
    data["gamma"] = kwargs[:alpha]

    # use beta if defined in setting or use 2*alpha
    # if haskey(kwargs, :beta)
    #     if kwargs[:beta] != 0
    #         data["beta"] = kwargs[:beta]
    #     end
    # else
    #     data["beta"] = 2*kwargs[:alpha]
    # end
    # # use gamma if defined in setting or use alpha
    # if haskey(kwargs, :gamma)
    #     if kwargs[:gamma] != 0
    #         data["gamma"] = kwargs[:gamma]
    #     end
    # else
    #     data["gamma"] = kwargs[:alpha]
    # end


end

## build method for Distributed PowerModel using ADMM algorithm
function build_dopf_app(pm::AbstractPowerModel)

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

    objective_min_fuel_and_consensus!(pm, objective_app!)
end

## method to set the APP algorithm objective
function objective_app!(pm::AbstractPowerModel)

    ## APP parameters
    alpha = pm.data["alpha"]
    beta  = pm.data["beta"]
    gamma = pm.data["gamma"]

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
                        v_neighbor = primal_variable[area][comp][ids][vstring]
                        v_local = primal_variable[area_id][comp][ids][vstring]
                        v_dual = dual_variable[area][comp][ids][vstring]

                        objective += beta/2 *  (v - v_local)^2 + gamma * v * (v_local - v_neighbor) + v * v_dual
                    end
                end
            end
        end
    end

    JuMP.@objective(pm.model, Min,  objective)
end

## method to update the dual variable value
function update_app!(data::Dict{String, <:Any})

    ## APP parameters
    alpha = data["alpha"]

    ## data
    area_id = string(get_area_id(data))
    primal_variable = data["shared_primal"]
    dual_variable = data["shared_dual"]

    ## update dual variable
    for area in keys(dual_variable)
        for comp in keys(dual_variable[area])
            for ids in keys(dual_variable[area][comp])
                for vstring in keys(dual_variable[area][comp][ids])

                    v_neighbor = primal_variable[area][comp][ids][vstring]
                    v_local = primal_variable[area_id][comp][ids][vstring]
                    v_dual = dual_variable[area][comp][ids][vstring]

                    data["shared_dual"][area][comp][ids][vstring] = v_dual  + alpha * (v_local - v_neighbor)
                end
            end
        end
    end

end
