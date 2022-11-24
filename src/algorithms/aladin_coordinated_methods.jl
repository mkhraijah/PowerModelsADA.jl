###############################################################################
#            Build methods for ADMM algorithm with coordinator                #
###############################################################################



# """
# ALADIN algorithm module contians build and update methods
# """

module aladin_coordinated_methods

using JuMP
using SparseArrays
using LinearAlgebra
using ..PowerModelsADA

"inilitlize the ALADIN algorithm local area"
function initialize_method_local(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    area_id = get_area_id(data)
    areas_id = get_areas_id(data)
    deleteat!(areas_id, areas_id .== area_id) # remove the same area from the list of areas_id

    # initiate primal and dual shared variables

    data["shared_variable"] = initialize_shared_variable(data, model_type, area_id, 0, "shared_variable", "flat")

    data["shared_dual_variable"] = Dict{String, Dict{String, Any}}("0" => initialize_shared_variable(data, model_type, area_id, areas_id, "shared_dual_variable", "flat"))

    data["received_dual_variable"] = Dict{String, Dict{String, Any}}("0" => initialize_shared_variable(data, model_type, area_id, areas_id, "received_dual_variable", "flat"))

    data["shared_sensitivities"] = Dict{String, Dict{String, Any}}("0" => Dict("g"=>Dict{String, Any}(), "C"=>Dict{String, Any}(), "B"=> Dict{String, Any}() ))

    data["local_solution"] = initialize_all_variable(data, model_type)

    data["received_delta"] = Dict{String, Dict{String, Any}}("0" => initialize_all_variable(data, model_type, "received_delta", "zeros"))

    # initiate algorithm settings
    initialize_dopf!(data, model_type; kwargs...)
    
    # initiate ALADIN parameters
    p = get(kwargs, :p, 1000)
    r_p = get(kwargs, :r_p, 1.3)
    p_upper = get(kwargs, :p_upper, 1e6)

    a1 = get(kwargs, :a1, 1)
    a2 = get(kwargs, :a2, 1)
    a3 = get(kwargs, :a3, 1)
    
    q_gamma = get(kwargs, :q_gamma, 0)
    sigma = get(kwargs, :sigma, NaN)

    data["parameter"] = Dict("p" => p, "r_p"=> r_p, "p_upper"=>p_upper, "a1"=>a1, "a2"=>a2, "a3"=> a3, "q_gamma"=>q_gamma, "sigma"=>sigma)

end


"inilitlize the ALADIN algorithm coordinator"
function initialize_method_coordinator(data::Dict{String, <:Any}, model_type::DataType; kwargs...)

    data_system = data
    data = deepcopy(data_system)
    data = decompose_coordinator(data)
    areas_id = get_areas_id(data)
    # initiate primal and dual shared variables

    data["received_sensitivities"] = Dict{String,Any}([string(area) => Dict{String, Any}() for area in areas_id])
    data["received_variable"] = initialize_shared_variable(data, model_type, 0 ,areas_id, "received_variable", "flat")
    data["shared_dual_variable"] = Dict{String,Any}()
    data["received_dual_variable"] = Dict{String,Any}()
    data["shared_delta"] = Dict{String,Any}()
    for i in areas_id
        data_area = decompose_system(data_system, i)
        areas = deepcopy(areas_id)
        deleteat!(areas, areas .== i)
        data["shared_dual_variable"][string(i)] = initialize_shared_variable(data, model_type, i ,areas, "shared_dual_variable", "flat")
        data["received_dual_variable"][string(i)] = initialize_shared_variable(data, model_type, i ,areas, "received_dual_variable", "flat")
        data["shared_delta"][string(i)] = initialize_all_variable(data_area, model_type, "shared_delta")
    end

    # initiate distributed algorithm parameters
    initialize_dopf!(data, model_type; kwargs...)

    mu = get(kwargs, :mu, 1000)
    r_mu = get(kwargs, :r_mu, 2)
    mu_upper = get(kwargs, :mu_upper, 2e6)

    a1 = get(kwargs, :a1, 1)
    a2 = get(kwargs, :a2, 1)
    a3 = get(kwargs, :a3, 1)
    
    q_gamma = get(kwargs, :q_gamma, 0)
    sigma = get(kwargs, :sigma, NaN)

    data["parameter"] = Dict("mu" => mu, "r_mu" => r_mu, "mu_upper" => mu_upper, "a1" => a1, "a2" => a2, "a3" => a3, "q_gamma" => q_gamma, "sigma" => sigma)

    return data
end

function calc_mismatch_aladin!(data::Dict{String, <:Any}; p::Int64=2 )
    
    area_id = string(get_area_id(data)) 
    mismatch_method = data["option"]["mismatch_method"]
    shared_variable_local = data["shared_variable"]["0"]
    shared_variable_solution = data["local_solution"]

    mismatch = Dict{String, Any}([
        variable => Dict{String, Any}([
            idx => shared_variable_local[variable][idx] - shared_variable_solution[variable][idx]
        for idx in keys(shared_variable_local[variable])]) 
    for variable in keys(shared_variable_local)])
    

    if mismatch_method == "norm"
        mismatch[area_id] = LinearAlgebra.norm([value for variable in keys(mismatch) for (idx,value) in mismatch[variable]], p)
    elseif mismatch_method == "max" || mismatch_method == "maximum"
        mismatch[area_id] = LinearAlgebra.maximum([value for variable in keys(mismatch) for (idx,value) in mismatch[variable]])
    end

    data["mismatch"] = mismatch
end

"update the ALADIN algorithm coordinator data before each iteration"
function update_method_local(data::Dict{String, <:Any})
   
    # parameters
    p = data["parameter"]["p"]
    r_p = data["parameter"]["r_p"]
    p_upper = data["parameter"]["p_upper"]
    a1 = data["parameter"]["a1"]
    a2 = data["parameter"]["a2"]
    a3 = data["parameter"]["a3"]

    ## data
    dual_variable_local = data["shared_dual_variable"]["0"]
    dual_variable_central = data["received_dual_variable"]["0"]
    delta = data["received_delta"]["0"]

    solution = data["solution"]
    local_solution = data["local_solution"]

    ## update dual variable
    for area in keys(dual_variable_local)
        for variable in keys(dual_variable_local[area])
            for idx in keys(dual_variable_local[area][variable])
                dual_central = dual_variable_central[area][variable][idx]
                dual_local = dual_variable_local[area][variable][idx]

                dual_variable_local[area][variable][idx]= dual_local + a3 * (dual_central - dual_local)
            end
        end
    end

    ## update solution (corresponding to update all variables)
    for variable in keys(local_solution)
        for idx in keys(local_solution[variable])
            v_sol = solution[variable][idx]
            v_local = local_solution[variable][idx]
            v_delta = delta[variable][idx]
            local_solution[variable][idx] = v_local + a1 * (v_sol - v_local) + a2 * v_delta
        end
    end


    # update parameters
    if p < p_upper
        data["parameter"]["p"] = r_p * p
    end

    PowerModelsADA.save_solution!(data)
    calc_mismatch_aladin!(data)
    PowerModelsADA.update_flag_convergence!(data)
    PowerModelsADA.update_iteration!(data)

end


"build PowerModel for ALADIN algorithm local area"
function build_method_local(pm::PowerModelsADA.AbstractPowerModel)

    # define variables
    variable_opf(pm)

    # define constraints
    constraint_opf(pm)

    # define objective function
    objective_min_fuel_and_consensus!(pm, objective_aladin_local)
end


"ALADIN algorithm objective coordinator"
function objective_aladin_local(pm::AbstractPowerModel)

    ## ALADIN parameters
    p = pm.data["parameter"]["p"]
    q_gamma = pm.data["parameter"]["q_gamma"]
    sigma = pm.data["parameter"]["sigma"]

    ## data
    area_id = get_area_id(pm)
    dual_variable = pm.data["shared_dual_variable"]["0"]
    local_solution = pm.data["local_solution"]

    ## objective function
    if haskey(pm.data["solution"], "qg")
        objective = q_gamma*(sum(PowerModelsADA._var(pm, "qg", string(idx)) for idx in ids(pm, :gen)))^2
    else
        objective = 0
    end


    for area in keys(dual_variable)
        for variable in keys(dual_variable[area])
            for idx in keys(dual_variable[area][variable])
                v = PowerModelsADA._var(pm, variable, idx)
                v_dual = dual_variable[area][variable][idx]
                if area_id < parse(Int64,area) 
                    objective += v * v_dual
                else
                    objective -= v * v_dual
                end
            end
        end
    end

    for variable in keys(local_solution)
        for idx in keys(local_solution[variable])
            v = PowerModelsADA._var(pm, variable, idx)
            v_local = local_solution[variable][idx]
            objective += sigma[variable] * p / 2 * (v - v_local)^2
        end
    end

    return objective
end

function update_sensitivities!(pm::AbstractPowerModel, solution::Dict{String, <:Any})

    solution["shared_sensitivities"] = pm.data["shared_sensitivities"]
    solution["shared_sensitivities"]["0"]["g"] = compute_cost_gradient(pm)
    solution["shared_sensitivities"]["0"]["C"] = compute_active_constraint_jacobian(pm, true)
    solution["shared_sensitivities"]["0"]["B"] = compute_optimal_hessian(pm)

end

post_processors_local = [update_solution!, update_shared_variable!, update_sensitivities!]

function compute_cost_gradient(pm::AbstractPowerModel)
    g = initialize_all_variable(pm.data, typeof(pm),"zeros")
    for (i,gen) in pm.data["gen"]
        pg = PowerModelsADA._var(pm, "pg", i)
        g["pg"][i] = 2*gen["cost"][1]*value(pg) + gen["cost"][2]
    end
    return g 
end

function compute_active_constraint_jacobian(pm::AbstractPowerModel, delete_dependent_constraints::Bool)
    # obtain the Jacobian matrix
    C_matrix = active_constraint_jacobian(pm, true)

    # arrange Jacobian matrix for each variable
    variable_dict = initialize_all_variable(pm.data, typeof(pm),"zeros")
    C = Dict([variable => Dict([idx=> Vector() for idx in keys(variable_dict[variable])]) for variable in keys(variable_dict)])
    for variable in keys(C)
        for idx in keys(C[variable])
            v = PowerModelsADA._var(pm, variable, idx)
            col = JuMP.index(v).value
            C[variable][idx] = C_matrix[:,col]
        end
    end
    return C
end


function active_constraint_jacobian(pm::AbstractPowerModel, delete_dependent_constraints::Bool)
    
    model = pm.model
    x = JuMP.all_variables(model)
    x_optimal = JuMP.value.(x)
    C = zeros(0,length(x))
    if !(isempty(JuMP.all_nonlinear_constraints(model)))
        d = JuMP.NLPEvaluator(model)
        JuMP.MOI.initialize(d, [:Jac])
        jacobian_sparsity = JuMP.MOI.jacobian_structure(d)
        I = [i for (i, _) in jacobian_sparsity]
        J = [j for (_, j) in jacobian_sparsity]
        V = zeros(length(jacobian_sparsity))
        nc = JuMP.num_nonlinear_constraints(model)
        JuMP.MOI.eval_constraint_jacobian(d, V, x_optimal)
        C = SparseArrays.sparse(I, J, V, nc, length(x))
    end
    
    #Now add in gradients for linear and quadratic constraints, which were not included in NLPEvaluator
    cref_types = JuMP.list_of_constraint_types(model)
    for t in cref_types
        if t[1] == JuMP.AffExpr   #Add gradients of all active linear constraints 
            crefs = JuMP.all_constraints(model, t[1], t[2])
            for cref in crefs
                if abs(JuMP.dual(cref)) > 1e-3 || t[2] == JuMP.MOI.EqualTo{Float64}  #only add gradient of active constraints 
                    cgrad = zeros(1, length(x))
                    for k = 1:length(x)
                        cgrad[1,k] = JuMP.normalized_coefficient(cref, x[k])
                    end
                    C = vcat(C, cgrad)
                end
            end
        elseif t[1] == QuadExpr  #Add gradients of all active quadratic constraints
            vmap = Dict(x[i] => i for i in 1:length(x))
            crefs = JuMP.all_constraints(model, t[1], t[2])
            for cref in crefs 
                if abs(dual(cref)) > 1e-3 || t[2] == JuMP.MOI.EqualTo{Float64}
                    cgrad = zeros(1,length(x))
                    func_set = JuMP.constraint_object(cref)
                    for (vars, coef) in func_set.func.terms
                        if vars.a != vars.b
                            cgrad[1,vmap[vars.a]] += x_optimal[vmap[vars.b]]*coef
                            cgrad[1,vmap[vars.b]] += x_optimal[vmap[vars.a]]*coef 
                        else
                            cgrad[1,vmap[vars.a]] += 2*x_optimal[vmap[vars.a]]*coef
                        end
                    end
                    for k = 1:length(x)
                        cgrad[1,k] += JuMP.coefficient(func_set.func, x[k])  #add coefficient of affine terms 
                    end
                    C = vcat(C, cgrad)
                end
            end
        end
    end
    #Add gradient of active variable bound constraints
    for k = 1:length(x)
        if JuMP.has_lower_bound(x[k])
            if abs(JuMP.value(x[k]) - JuMP.lower_bound(x[k])) < 1e-6
                cgrad = zeros(1, length(x))
                cgrad[1,k] = 1
                C = vcat(C, cgrad)
            end
        end
        if JuMP.has_upper_bound(x[k])
            if abs(JuMP.value(x[k]) - JuMP.upper_bound(x[k])) < 1e-6
                cgrad = zeros(1, length(x))
                cgrad[1,k] = 1
                C = vcat(C, cgrad)
            end
        end
        if JuMP.is_fixed(x[k])
            cgrad = zeros(1, length(x))
            cgrad[1,k] = 1
            C = vcat(C, cgrad)
        end
    end

    # Delete linearly dependent rows until constraints are all linearly independent (Ipopt doesn't do this, although Gurobi would)
    if (delete_dependent_constraints)
        C = delete_dependent_row(C)
    end

    return Matrix(C)
end

function delete_dependent_row(mat)
    while(size(mat)[1] > LinearAlgebra.rank(mat))
        x = findfirst(i -> rank(mat[1:end .!= i, :]) == rank(mat), 1:size(mat)[1])
        mat =mat[1:end .!= x, :]
    end
    return mat
end

function compute_optimal_hessian(pm::AbstractPowerModel)
    # obtain the Hessian matrix
    B_matrix = optimal_hessian(pm)

    # arrange Hessian matrix for each variable
    variable_dict = initialize_all_variable(pm.data, typeof(pm),"zeros")
    B = Dict([variable1 => Dict([idx1=> Dict(
        [variable2 => Dict([idx2=> 0.0 for idx2 in keys(variable_dict[variable2])]) for variable2 in keys(variable_dict)]
    ) for idx1 in keys(variable_dict[variable1])]) for variable1 in keys(variable_dict)])
    for variable1 in keys(B)
        for idx1 in keys(B[variable1])
            v1 = PowerModelsADA._var(pm, variable1, idx1)
            for variable2 in keys(B)
                for idx2 in keys(B[variable2])
                    v2 = PowerModelsADA._var(pm, variable2, idx2)
                    row = JuMP.index(v1).value
                    col = JuMP.index(v2).value
                    B[variable1][idx1][variable2][idx2] = B_matrix[row,col]
                end
            end
        end
    end
    return B
end

function optimal_hessian(pm::AbstractPowerModel)

    pm_temp = deepcopy(pm)
    model = pm_temp.model
    x = all_variables(model)
    n = num_variables(model)
    x_optimal = value.(x)
    y_optimal = -dual.(all_nonlinear_constraints(model))
    crefs = Dict{ConstraintRef,Float64}()
    for (F, S) in list_of_constraint_types(model)
        for cref in all_constraints(model, F, S)
            crefs[cref] = -dual(cref)
        end
    end

    #change the objective to generator cost only
    PowerModelsADA._PM.objective_min_fuel_and_flow_cost(pm_temp)

    # compute hessian matrix
    H = spzeros(length(x),length(x))
    if !(isempty(all_nonlinear_constraints(model)))
        d = NLPEvaluator(model)
        MOI.initialize(d, [:Hess])
        hessian_sparsity = MOI.hessian_lagrangian_structure(d)
        I = [i for (i, _) in hessian_sparsity]
        J = [j for (_, j) in hessian_sparsity]
        V = zeros(length(hessian_sparsity))
        MOI.eval_hessian_lagrangian(d, V, x_optimal, 1.0, y_optimal)
        H = SparseArrays.sparse(I, J, V, n, n)
    end
    vmap = Dict(x[i] => i for i in 1:n)
    for (F, S) in list_of_constraint_types(model)
        for cref in all_constraints(model, F, S)
            add_to_hessian(H, constraint_object(cref).func, crefs[cref], vmap)
        end
    end
    add_to_hessian(H, objective_function(model), 1.0, vmap)
    B = fill_off_diagonal(H)

    if findmin(eigvals(Matrix(B)))[1] < 1e-4
        F = eigen(Matrix(B));
        d = zeros(length(F.values));
        for idx = 1:length(F.values)
            if F.values[idx] < -1e-4
                d[idx] = abs(F.values[idx])
                d[idx] = -1*(F.values[idx])
            elseif -1e-4 <= F.values[idx] <= 1e-4
                d[idx] = 1e-4 
            else
                d[idx] = F.values[idx]
            end
        end

        D = LinearAlgebra.Diagonal(d)

        B = F.vectors*abs.(D)*transpose(F.vectors)
    end
    return B 

end

add_to_hessian(H, f::Any, μ, vmap) = nothing
function add_to_hessian(H, f::QuadExpr, μ, vmap)
    for (vars, coef) in f.terms
        if vars.a != vars.b 
            H[vmap[vars.a], vmap[vars.b]] += μ * coef
        else
            H[vmap[vars.a], vmap[vars.b]] += 2 * μ * coef
        end
    end
end

#helper function to fill in missing symmetric elements in sparse Hessian
function fill_off_diagonal(H)
    ret = H + H'
    row_vals = SparseArrays.rowvals(ret)
    non_zeros = SparseArrays.nonzeros(ret)
    for col in 1:size(ret, 2)
        for i in SparseArrays.nzrange(ret, col)
            if col == row_vals[i]
                non_zeros[i] /= 2
            end
        end
    end
    return ret
end

"build PowerModel for ALADIN algorithm coordinator"
function solve_coordinator!(data, optimizer)
    
    mu = data["parameter"]["mu"]

    shared_variable = data["received_variable"]
    dual_variable = data["received_dual_variable"]
    sensitivities = data["received_sensitivities"]
    delta = data["shared_delta"]
    qp_dual_variable = data["shared_dual_variable"]

    model = Model(optimizer)

    # define variables
    x = Dict{String, Any}(area => Dict{String, Any}(variable => Dict{String, Any}(idx => JuMP.@variable(model, base_name=string("x_", area, "_", variable, "_", idx)) for idx in keys(delta[area][variable])) for variable in keys(delta[area])) for area in keys(delta))

    s = Dict{String, Any}(area1 => Dict{String, Any}(area2=> Dict{String, Any}(variable => Dict{String, Any}(idx => JuMP.@variable(model, base_name=string("s_", variable, "_", idx)) for idx in keys(dual_variable[area1][area2][variable])) for variable in keys(dual_variable[area1][area2])) for area2 in keys(dual_variable[area1]) if area1<area2) for area1 in keys(dual_variable))


    # define objective function
    qp_objective = JuMP.GenericQuadExpr{Float64, JuMP.VariableRef}()
    for area1 in keys(dual_variable)
        for area2 in keys(dual_variable[area1])
            if area1 < area2
                for variable in keys(dual_variable[area1][area2])
                    for (idx,val) in dual_variable[area1][area2][variable]
                        qp_objective += val*s[area1][area2][variable][idx] + mu/2*(s[area1][area2][variable][idx])^2
                    end
                end
            end
        end
    end

    for area in keys(x)
        for variable in keys(x[area])
            for idx in keys(x[area][variable])
                qp_objective += sum(0.5*x[area][variable][idx]*x[area][variable2][idx2]*sensitivities[area]["B"][variable][idx][variable2][idx2] for variable2 in keys(x[area]) for idx2 in keys(x[area][variable2]))+ (sensitivities[area]["g"][variable][idx]*x[area][variable][idx])
            end
        end
    end

    @objective(model, Min, qp_objective)

    # define constraints
    constraint_ref = Dict{String,Any}(area1 => Dict{String,Any}(area2 => Dict{String,Any}(variable => Dict{String, JuMP.ConstraintRef}([idx => @constraint(model, x[area1][variable][idx] + shared_variable[area1][variable][idx] - x[area2][variable][idx] - shared_variable[area2][variable][idx] == s[area1][area2][variable][idx]) for idx in keys(dual_variable[area1][area2][variable])]) for variable in keys(dual_variable[area1][area2]) ) for area2 in keys(dual_variable[area1]) if area1<area2) for area1 in keys(dual_variable))


    for area in keys(x)
        n_constraints = size(first(first(data["received_sensitivities"][area]["C"])[2])[2])[1]
        for n in 1:n_constraints
            @constraint(model, sum(sensitivities[area]["C"][variable][idx][n]*x[area][variable][idx] for variable in keys(x[area]) for idx in keys(x[area][variable])) == 0)
        end
    end


    
    optimize!(model)

    for area in keys(delta)
        for variable in keys(delta[area])
            for idx in keys(delta[area][variable])
                delta[area][variable][idx] = value(x[area][variable][idx])
            end
        end
    end

    for area1 in keys(qp_dual_variable)
        for area2 in keys(qp_dual_variable[area1])
            if area1 < area2
                for variable in keys(qp_dual_variable[area1][area2])
                    for idx in keys(qp_dual_variable[area1][area2][variable])
                        qp_dual_variable[area1][area2][variable][idx] = -dual(constraint_ref[area1][area2][variable][idx])
                        qp_dual_variable[area2][area1][variable][idx] = -dual(constraint_ref[area1][area2][variable][idx])
                    end
                end
            end
        end
    end

end


"update the ADMM algorithm coordinator data after each iteration"
function update_method_coordinator(data::Dict{String, <:Any}) 

    if data["parameter"]["mu"] < data["parameter"]["mu_upper"]
        data["parameter"]["mu"] = data["parameter"]["r_mu"]*data["parameter"]["mu"]
    end

    save_solution!(data)
    update_iteration!(data)

end

end




"""
    solve_dopf_aladin_coordinated(data::Dict{String, <:Any}, model_type::DataType, optimizer; tol::Float64=1e-4, 
    max_iteration::Int64=1000, print_level = true, p::Real=1000, mu::Real=1000, p_upper::Real=1e6, mu_upper::Real=2e6, r_p::Real=1.5, mu_p::Real=2, a1::Real=1, a2::Real=1, #     a3::Real=1, q_gamma::Real=0, sigma::Dict{String,Real}=Dict())

Solve the distributed OPF problem using ADMM algorithm with central coordinator.

# Arguments:
- data::Dict{String, <:Any} : dictionary contains case in PowerModel format
- model_type::DataType : power flow formulation (PowerModel type)
- optimizer : optimizer JuMP initiation object
- mismatch_method::String="norm" : mismatch calculation method (norm, max)
- tol::Float64=1e-4 : mismatch tolerance
- max_iteration::Int64=1000 : maximum number of iteration
- print_level::Int64=1 : print mismatch after each iteration and result summary 
- p::Real=1000 : algorithm parameters
"""
function solve_dopf_aladin_coordinated(data::Union{Dict{String, <:Any}, String}, model_type::DataType, optimizer; mismatch_method::String="norm", tol::Float64=1e-4, max_iteration::Int64=1000, print_level::Int64=1, p::Real=1000, mu::Real=1000, p_upper::Real=1e6, mu_upper::Real=2e6, r_p::Real=1.5, r_mu::Real=2, a1::Real=1, a2::Real=1, a3::Real=1, q_gamma::Real=0, sigma::Dict{String, <:Any}=Dict())


    # obtain and standardize case data
    if isa(data, String)
        data = parse_file(data)
    end
    _PM.standardize_cost_terms!(data, order=2)

    # arrange sigma parameter
    if isempty(sigma)
        sigma = Dict{String,Any}("w"=> 20, "wr"=>5, "wi"=>5 ,"vi"=> 10, "vr"=> 10, "va" => 10, "vm" => 5, "pf" => 1, "pt" => 1, "qf" => 1, "qt" => 1, "pg" => 1, "qg" => 1)
    end

    # obtain and arrange areas id
    arrange_areas_id!(data)
    areas_id = get_areas_id(data)

    # decompose the system into subsystems
    data_area = Dict{Int64, Any}()
    for i in areas_id
        data_area[i] = decompose_system(data, i)
    end

    # initilize distributed power model parameters
    data_coordinator = aladin_coordinated_methods.initialize_method_coordinator(data, model_type; mismatch_method=mismatch_method, max_iteration=max_iteration, tol=tol, a1=a1, a2=a2, a3=a3, mu=mu, r_mu=r_mu, mu_upper=mu_upper, q_gamma=q_gamma, sigma=sigma)
    for i in areas_id
        aladin_coordinated_methods.initialize_method_local(data_area[i], model_type; mismatch_method=mismatch_method, max_iteration=max_iteration, tol=tol, p=p, a1=a1, a2=a2, a3=a3, r_p=r_p, p_upper=p_upper, q_gamma=q_gamma, sigma=sigma)
    end

    ## initialaize the algorithms global counters
    iteration = 0
    flag_convergence = false

    ## start iteration
    while iteration < max_iteration && !flag_convergence

        # solve local area problems in parallel
        info1 = @capture_out begin
            Threads.@threads for i in areas_id
                result = solve_model(data_area[i], model_type, optimizer, aladin_coordinated_methods.build_method_local, solution_processors=aladin_coordinated_methods.post_processors_local)
                update_data!(data_area[i], result["solution"])
            end
        end

        # share solution of local areas with the coordinator
        for i in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_area[i], 0, serialize = false)
            receive_shared_data!(data_coordinator, shared_data, i)
        end

        # solve coordinator problem 
        info2 = @capture_out begin
            aladin_coordinated_methods.solve_coordinator!(data_coordinator, optimizer)
        end

        # share coordinator solution with local areas
        for i in areas_id # sender subsystem
            shared_data = prepare_shared_data(data_coordinator, i, serialize = false)
            receive_shared_data!(data_area[i], shared_data, 0)
        end

        # update local areas and coordinator problems after
        aladin_coordinated_methods.update_method_coordinator(data_coordinator)
        for i in areas_id
            aladin_coordinated_methods.update_method_local(data_area[i])
        end

        # check global convergence and update iteration counters
        flag_convergence = update_global_flag_convergence(data_area)
        iteration += 1

        # print solution
        print_iteration(data_area, print_level, [info1; info2])

    end

    data_area[0] = data_coordinator
    print_convergence(data_area, print_level)

    return data_area
end

# export the algorithm methods module and call method
export aladin_coordinated_methods, solve_dopf_aladin_coordinated