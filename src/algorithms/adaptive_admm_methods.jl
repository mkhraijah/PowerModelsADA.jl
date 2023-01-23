###############################################################################
#                 Build methods for adaptive ADMM algorithm                   #
###############################################################################

"""
adaptive ADMM algorithm module contains build and update methods
"""
module adaptive_admm_methods
using ..PowerModelsADA
using LinearAlgebra
"solve distributed OPF using adaptive ADMM algorithm"
function solve_method(data, model_type::DataType, optimizer; kwargs...)
    solve_dopf(data, model_type, optimizer, adaptive_admm_methods; kwargs...)
end

"initialize the adaptive ADMM algorithm"
function initialize_method(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    area_id =get_area_id(data)
    areas_id = get_areas_id(data)
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id
    initialization_method = get(kwargs, :initialization_method, "flat")

    # primal and dual shared variables
    data["shared_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "shared_variable", initialization_method)

    data["received_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "received_variable", initialization_method)

    data["dual_variable"] = initialize_shared_variable(data, model_type, area_id, areas_id, "dual_variable", initialization_method)

    data["dual_residual"] = Dict{String, Any}()
    # distributed algorithm settings
    initialize_dopf!(data, model_type; kwargs...)

    if haskey(data, "previous_solution")
        for str in unique(["shared_variable", "received_variable",keys(data["previous_solution"])...])
            if !haskey(data["previous_solution"], str)
                data["previous_solution"][str]= Vector{Dict}()
            end
        end
    else
        data["previous_solution"]= Dict{String, Any}([str=> Vector{Dict}() for str in ["shared_variable", "received_variable"]])
    end

    # adaptive ADMM parameters
    alpha = Float64(get(kwargs, :alpha, 1000))
    mu_inc = Float64(get(kwargs, :mu_inc, 2))
    mu_dec = Float64(get(kwargs, :mu_dec, 2))
    eta_inc = Float64(get(kwargs, :eta_inc, 0.2))
    eta_dec = Float64(get(kwargs, :eta_dec, 0.2))

    data["parameter"] = Dict("alpha"=> alpha,"mu_inc"=> mu_inc, "mu_dec"=> mu_dec, "eta_inc"=> eta_inc, "eta_dec"=>eta_dec)
    data["alpha"] = initialize_shared_variable(data, model_type, area_id, areas_id, "parameter", "constant", alpha)
end

"build PowerModel object for the adaptive ADMM algorithm"
function build_method(pm::AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)
  
    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_adaptive_admm)
end

"adaptive ADMM algorithm objective function"
function objective_adaptive_admm(pm::AbstractPowerModel)

    # parameters
    alphas = pm.data["alpha"]

    # data
    shared_variable_local = pm.data["shared_variable"]
    shared_variable_received = pm.data["received_variable"]
    dual_variable = pm.data["dual_variable"]

    ##objective function
    objective = 0
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PowerModelsADA._var(pm, variable, idx)
                v_central = (shared_variable_local[area][variable][idx] + shared_variable_received[area][variable][idx])/2
                v_dual = dual_variable[area][variable][idx]
                alpha = alphas[area][variable][idx]
                objective += alpha/2 * (v - v_central)^2 + v_dual * (v - v_central)
            end
        end
    end

    return objective
end

"update the adaptive ADMM algorithm data after each iteration"
function update_method(data::Dict{String, <:Any})
    
    # parameters
    alphas = data["alpha"]
    mu_inc = data["parameter"]["mu_inc"]
    mu_dec = data["parameter"]["mu_dec"]
    eta_inc = data["parameter"]["eta_inc"]
    eta_dec = data["parameter"]["eta_dec"]

    # data
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]
    dual_variable = deepcopy(data["dual_variable"])

    # update dual variable
    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v_primal = shared_variable_local[area][variable][idx]
                v_central = (shared_variable_local[area][variable][idx] + shared_variable_received[area][variable][idx])/2
                v_dual = dual_variable[area][variable][idx]
                alpha = alphas[area][variable][idx]
                data["dual_variable"][area][variable][idx]= v_dual + alpha * (v_primal - v_central)
            end
        end
    end

    calc_dual_residual_adaptive!(data)
    calc_mismatch!(data, central = true)

    ## update adaptive ADMM parameters
    for area in keys(data["alpha"])
        for variable in keys(data["alpha"][area])
            for idx in keys(data["alpha"][area][variable])
                if data["mismatch"][area][variable][idx] > mu_inc * data["dual_residual"][area][variable][idx]
                    data["alpha"][area][variable][idx] = data["alpha"][area][variable][idx] * ( 1 + eta_inc)
                elseif data["dual_residual"][area][variable][idx] > mu_dec * data["mismatch"][area][variable][idx]
                    data["alpha"][area][variable][idx] = data["alpha"][area][variable][idx] / ( 1 + eta_dec)
                end
            end
        end
    end

    update_flag_convergence!(data)
    save_solution!(data)
    update_iteration!(data)
end


"""
    calc_dual_residual!(data::Dict{String, <:Any}; central::Bool=false)
calculate the dual redidual as seen by the area. Set central=true if the algorithm uses the optimality condition of a central coordinator.
"""
function calc_dual_residual_adaptive!(data::Dict{String, <:Any}; central::Bool=false)
    area_id = string(get_area_id(data))
    alpha = data["alpha"]
    shared_variable_local = data["shared_variable"]
    shared_variable_received = data["received_variable"]

    if data["counter"]["iteration"] == 1
        dual_dual_residual = Dict{String, Any}([
            area => Dict{String, Any}([
                variable => Dict{String, Any}([
                    idx => central ? -alpha[area][variable][idx]* (shared_variable_local[area][variable][idx]+shared_variable_received[area][variable][idx])/2 : -alpha[area][variable][idx]* shared_variable_local[area][variable][idx]
                for idx in keys(shared_variable_local[area][variable])])
            for variable in keys(shared_variable_local[area])])
        for area in keys(shared_variable_local)])
    else
        previous_shared_variable_local = data["previous_solution"]["shared_variable"][end]
        previous_shared_variable_received = data["previous_solution"]["received_variable"][end]
        dual_dual_residual = Dict{String, Any}([
            area => Dict{String, Any}([
                variable => Dict{String, Any}([
                    idx => central ? -alpha[area][variable][idx] * ((shared_variable_local[area][variable][idx]+shared_variable_received[area][variable][idx])/2 - (previous_shared_variable_local[area][variable][idx] +previous_shared_variable_received[area][variable][idx] )/2) : -alpha[area][variable][idx] * (shared_variable_local[area][variable][idx] - previous_shared_variable_local[area][variable][idx])
                for idx in keys(shared_variable_local[area][variable])])
            for variable in keys(shared_variable_local[area])])
        for area in keys(shared_variable_local) ])
    end

    dual_dual_residual[area_id] = LinearAlgebra.norm([value for area in keys(dual_dual_residual) if area != area_id for variable in keys(dual_dual_residual[area]) for (idx,value) in dual_dual_residual[area][variable]])
    data["dual_residual"] = dual_dual_residual

end

post_processors = [update_solution!, update_shared_variable!]
end

"""
    solve_dopf_adaptive_admm(data::Dict{String, <:Any}, model_type::DataType, optimizer; 
    mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, 
    print_level::Int64=1, print_optimizer_info::Bool=false, alpha::Real=1000)

Solve the distributed OPF problem using adaptive ADMM algorithm.

# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- print_level::Int64=1 : 0 - no print, 1 - print mismatch after each iteration and result summary, 2 - print optimizer output
- alpha::Real=1000 : algorithm parameter
"""
solve_dopf_adaptive_admm = adaptive_admm_methods.solve_method

# export the algorithm methods module and solve method
export adaptive_admm_methods, solve_dopf_adaptive_admm