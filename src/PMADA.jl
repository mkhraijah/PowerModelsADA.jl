
module PMADA

import JuMP
import PowerModels as _PM
import Serialization
import LinearAlgebra
import DelimitedFiles
import Clustering: kmeans

import PowerModels: AbstractPowerModel, parse_file, ids, ref, var, con, sol, nw_ids, nws, optimize_model!, nw_id_default, ismultinetwork, update_data!


# PowerModels types
powermodels = names(_PM)
powermodels = filter(x -> endswith(string(x), "PowerModel"), powermodels)
powermodels = filter(x -> !occursin("Abstract", string(x)), powermodels)
for type in powermodels
    @eval import PowerModels: $(type)
    @eval export $(type)
end

include("core/base.jl")
include("core/variables.jl")
include("core/objective.jl")
include("core/data.jl")
include("core/data_sharing.jl")
include("core/util.jl")

include("algorithms/admm_methods.jl")
include("algorithms/atc_methods.jl")
include("algorithms/app_methods.jl")

include("algorithms/coordinated_admm_methods.jl")
include("algorithms/coordinated_atc_methods.jl")

export
parse_file,
assign_area!,
solve_dopf_admm,
solve_dopf_app,
solve_dopf_atc,
solve_dopf_admm_coordinated,
solve_dopf_atc_coordinated,
compare_solution,
partition_system!,
calc_dist_gen_cost,
solve_local!,
calc_mismatch!,
update_flag_convergance!,
update_iteration!,
update_global_flag_convergance,
calc_global_mismatch,
arrange_areas_id!,
get_areas_id,
get_area_id,
get_local_bus,
get_neighbor_bus,
get_areas_bus,
get_shared_component,
initialize_variable_shared!,
initialize_primal_variable_shared!,
initialize_dual_variable_shared!,
initialize_variable_shared_coordinator!,
initialize_primal_variable_shared_coordinator!,
initialize_dual_variable_shared_coordinator!,
initialize_variable_shared_local!,
initialize_primal_variable_shared_local!,
initialize_dual_variable_shared_local!,
variable_shared_names,
objective_min_fuel_and_consensus!,
send_shared_data,
receive_shared_data!,
initialize_dopf_admm!,
build_dopf_admm,
objective_admm!,
update_admm!,
initialize_dopf_admm_local!,
initialize_dopf_admm_coordinator!,
build_dopf_admm_local,
build_dopf_admm_coordinator,
objective_admm_coordinated!,
update_admm_coordinated!,
initialize_dopf_atc!,
build_dopf_atc,
objective_atc!,
update_atc!,
initialize_dopf_atc_local!,
initialize_dopf_atc_coordinator!,
build_dopf_atc_local,
build_dopf_atc_coordinator,
objective_atc_coordinated!,
update_atc_coordinated!,
initialize_dopf_app!,
build_dopf_app,
objective_app!,
update_app!

# JuMP optimizer call
import JuMP: optimizer_with_attributes
export optimizer_with_attributes


export ids, ref, var, con, sol, nw_ids, nws, optimize_model!, nw_id_default, ismultinetwork, update_data!

end