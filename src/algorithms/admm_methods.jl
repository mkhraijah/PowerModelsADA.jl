###############################################################################
#                     Build methods for ADMM algorithm                        #
###############################################################################

"""
    solve_dopf_admm(data::Dict{String, <:Any}, model_type::Type, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, alpha::Real=1000)

Solve the distributed OPF problem using ADMM algorithm.
"""
function solve_dopf_admm(data::Dict{String, <:Any}, model_type::Type, optimizer; tol::Float64=1e-4, max_iteration::Int64=1000, verbose = true, alpha::Real=1000)

    solve_dopf(data, model_type, build_dopf_admm, update_admm!, optimizer, initialize_method=initialize_dopf_admm!, tol = tol, max_iteration=max_iteration, verbose=verbose, alpha=alpha)

end

"inilitlize the admm algorithm"
function initialize_dopf_admm!(data::Dict{String, <:Any}, model_type::Type; tol::Float64=1e-4, max_iteration::Int64=1000, kwargs)

    initialize_dopf!(data, model_type, tol=tol, max_iteration=max_iteration, kwargs=kwargs)
    data["alpha"] = kwargs[:alpha]

end

"build PowerModel using ADMM algorithm"
function build_dopf_admm(pm::AbstractPowerModel)

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

    objective_min_fuel_and_consensus!(pm, objective_admm!)
end

"ADMM algorithm objective"
function objective_admm!(pm::AbstractPowerModel)

    ## ADMM parameters
    alpha = pm.data["alpha"]

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

                        objective += alpha/2 * (v - v_central)^2 + v_dual * (v - v_central)
                    end
                end
            end
        end
    end

    JuMP.@objective(pm.model, Min,  objective)
end



"update the dual variable value"
function update_admm!(data::Dict{String, <:Any})

    ## ADMM parameters
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

                    v_primal = primal_variable[area_id][comp][ids][vstring]
                    v_central = (primal_variable[area_id][comp][ids][vstring] + primal_variable[area][comp][ids][vstring])/2
                    v_dual = dual_variable[area][comp][ids][vstring]

                    data["shared_dual"][area][comp][ids][vstring] = v_dual  + alpha * (v_primal - v_central)
                end
            end
        end
    end
end
