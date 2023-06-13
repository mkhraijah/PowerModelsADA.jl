###############################################################################
#       Methods for data wrangling and extracting information from data       #
###############################################################################


"assign area to the system data using a dictionary with (bus => area) Int pairs"
function assign_area!(data::Dict{String, <:Any}, partition::Dict)
    for i in keys(data["bus"])
        data["bus"][i]["area"] = partition[parse(Int64,i)]
    end
end

"assign area to the system data using a CVS file with buses and area ID"
function assign_area!(data::Dict{String, <:Any}, partition_path::String)
    partition_mat = DelimitedFiles.readdlm(partition_path, ',', Int, '\n')
    # explicit conversion from Matrix{Int64} to Vector{Pair{Int64, Int64}}
    partition = Vector{Pair{Int64, Int64}}([
        Pair{Int64, Int64}(row[1], row[2])
        for row in eachrow(partition_mat)
    ])
    assign_area!(data, partition)
end

"assign area to the system data using a vector with (bus => area) pairs"
function assign_area!(data::Dict{String, <:Any}, partition::Vector{Pair{Int64, Int64}})
    assign_area!(data, Dict(partition))
end

"assign area to the system data using a matrix with [bus, area] columnsor rows"
function assign_area!(data::Dict{String, <:Any}, partition::Array{Int64, 2})
    if size(partition)[2] != 2 && length(data["bus"]) != 2
        partition = partition'
        if size(partition)[2] != 2
            #through error
            error("Partitioning data does not contain correct area assignments")
        end
    end
    assign_area!(data, Dict(partition[i,1] => partition[i,2] for i in 1:size(partition)[1] ))
end

"""
    decompose_system(data::Dict{String, <:Any})

decompose a system into areas defined by bus area.
"""
function decompose_system(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    data_area = Dict{Int64, Any}([i => decompose_system(data, i) for i in areas_id])
    return data_area
end

"obtain an area decomposition with area ID"
function decompose_system(data::Dict{String, <:Any}, area_id::Int64)
    # identify local buses
    local_bus = get_local_bus(data, area_id)
    neighbor_bus = get_neighbor_bus(data, area_id)

    ## add virtual generators
    virtual_gen = add_virtual_gen(data, neighbor_bus, area_id)

    ## area data
    data_area = Dict{String,Any}()
    data_area["area"] = area_id
    data_area["name"]= "$(data["name"])_area_$area_id"
    data_area["source_version"] = data["source_version"]
    data_area["source_type"] = data["source_type"]
    data_area["baseMVA"] = data["baseMVA"]
    data_area["per_unit"] = data["per_unit"]
    data_area["bus"] = Dict{String,Any}([j => bus for (j,bus) in data["bus"] if bus["bus_i"] in [local_bus;neighbor_bus] && bus["bus_type"] != 4])
    data_area["branch"] = Dict{String,Any}([j => branch for (j,branch) in data["branch"] if (branch["f_bus"] in local_bus || branch["t_bus"] in local_bus) && branch["br_status"] == 1])
    data_area["gen"] = merge(Dict{String, Any}([i => gen for (i,gen) in data["gen"] if gen["gen_bus"] in local_bus && gen["gen_status"] == 1]), virtual_gen)
    data_area["shunt"] = Dict{String, Any}([i => shunt for (i,shunt) in data["shunt"] if shunt["shunt_bus"] in local_bus])
    data_area["load"] = Dict{String, Any}([i => load for (i,load) in data["load"] if load["load_bus"] in local_bus])
    data_area["storage"]= Dict{String, Any}([i => storage for (i,storage) in data["storage"] if gen["storage_bus"] in local_bus])
    data_area["switch"]=Dict{String, Any}([i => switch for (i,switch) in data["switch"] if gen["switch_bus"] in local_bus])
    data_area["dcline"]= Dict{String, Any}([i => dcline for (i,dcline) in data["dcline"] if dcline["f_bus"] in local_bus || dcline["t_bus"] in local_bus ] )

    return data_area
end

"obtain system coordinator data"
function decompose_coordinator(data::Dict{String, <:Any})
   
    areas_id = get_areas_id(data)

    # identifylocal buses
    boundary_bus = unique(reduce(vcat, [get_neighbor_bus(data, area_id) for area_id in areas_id]))

    ## area data
    data_coordinator = Dict{String,Any}()
    data_coordinator["area"] = 0
    data_coordinator["name"]= "$(data["name"])_coordinator"
    data_coordinator["source_version"] = data["source_version"]
    data_coordinator["source_type"] = data["source_type"]
    data_coordinator["baseMVA"] = data["baseMVA"]
    data_coordinator["per_unit"] = data["per_unit"]
    data_coordinator["bus"] = Dict{String,Any}([j => bus for (j,bus) in data["bus"] if bus["bus_i"] in boundary_bus])
    data_coordinator["branch"] = Dict{String,Any}([j => branch for (j,branch) in data["branch"] if branch["f_bus"] in boundary_bus && branch["t_bus"] in boundary_bus && data["bus"]["$(branch["f_bus"])"]["area"] != data["bus"]["$(branch["t_bus"])"]["area"] ])
    data_coordinator["gen"] = Dict{String,Any}()
    data_coordinator["shunt"] = Dict{String,Any}()
    data_coordinator["load"] = Dict{String,Any}()
    data_coordinator["storage"]= Dict{String,Any}()
    data_coordinator["switch"]= Dict{String,Any}()
    data_coordinator["dcline"]= Dict{String,Any}()
    return data_coordinator
end

"add virtual generators at the neighboring buses of an area"
function add_virtual_gen(data::Dict{String, <:Any}, neighbor_bus::Vector, area_id::Int64)
    max_gen_ind = maximum([parse(Int,i) for i in keys(data["gen"])])
    virtual_gen = Dict{String, Any}()
    cost_model = data["gen"][string(max_gen_ind)]["model"]
    max_flow = 10*sum(load["pd"] for (i,load) in data["load"])
    if cost_model == 1
        for i in neighbor_bus
            virtual_gen[string(i+max_gen_ind)] = Dict{String, Any}("ncost" => 2, "qc1max" => 0.0, "pg" => 0, "model" => cost_model, "shutdown" => 0.0, "startup" => 0.0, "qc2max" => 0.0, "ramp_agc" => 0.0, "qg" => 0.0, "gen_bus" => i, "pmax" => max_flow, "ramp_10" => 0.0, "vg" => 1.05, "mbase" => data["baseMVA"], "source_id" => Any["gen", max_gen_ind+i], "pc2" => 0.0, "index" => i+max_gen_ind, "cost" => [0.01; 0.0; 0.02; 0.0], "qmax" => max_flow, "gen_status" => 1, "qmin" => -max_flow, "qc1min" => 0.0, "qc2min" => 0.0, "pc1" => 0.0, "ramp_q" => 0.0, "ramp_30" => 0.0, "pmin" => -max_flow, "apf" => 0.0)
        end
    else
        for i in neighbor_bus
            virtual_gen[string(i+max_gen_ind)] = Dict{String, Any}("ncost" => 3, "qc1max" => 0.0, "pg" => 0, "model" => cost_model, "shutdown" => 0.0, "startup" => 0.0, "qc2max" => 0.0, "ramp_agc" => 0.0, "qg" => 0.0, "gen_bus" => i, "pmax" => max_flow, "ramp_10" => 0.0, "vg" => 1.05, "mbase" => data["baseMVA"], "source_id" => Any["gen", max_gen_ind+i], "pc2" => 0.0, "index" => i+max_gen_ind, "cost" => [0.0; 0.0; 0.0], "qmax" => max_flow, "gen_status" => 1, "qmin" => -max_flow, "qc1min" => 0.0, "qc2min" => 0.0, "pc1" => 0.0, "ramp_q" => 0.0, "ramp_30" => 0.0, "pmin" => -max_flow, "apf" => 0.0)
        end
    end
    return virtual_gen
end

# "add virtual bus at the tie-lines with other areas"
# function _add_virtual_bus!(data::Dict{String, <:Any}, neighbor_bus::Vector, area_id::Int)
#     max_bus_ind = maximum([parse(Int,i) for i in keys(data["bus"])])
#     vmax = first(data["bus"])[2]["vmax"]
#     vmin = first(data["bus"])[2]["vmin"]

#     virtual_bus = Dict{String, Any}()
#     for i in neighbor_bus
#         bus_area = data["bus"][string(i)]["area"]
#         base_kv = data["bus"][string(i)]["base_kv"]
#         common_lines = [idx for (idx,branch) in data["branch"] if (branch["f_bus"] == i && data["bus"][string(branch["t_bus"])]["area"] == area_id) || (branch["t_bus"] == i && data["bus"][string(branch["f_bus"])]["area"] == area_id) ]
#         for j in common_lines
#             bus_id = parse(Int64, j) + max_bus_ind
#             virtual_bus[string(bus_id)] = Dict{String, Any}("zone" => bus_area, "bus_i" => bus_id, "bus_type" => 1, "vmax" => vmax, "source_id" => Any["bus", bus_id], "area"=> bus_area, "vmin" => vmin, "index" => 0.0, "va" => 1.0, "vm" => 0.0,  "base_kv" => base_kv)
#         end
#     end
#     return virtual_bus
# end

"""
arrange area ID from 1 to number of areas. This step is necessary when having area number 0 and using central coordinator
"""
function arrange_areas_id!(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    new_areas_id = collect(1:length(areas_id))
    area_id_lookup = Dict(areas_id[i] => i for i in new_areas_id)
    for (i,bus) in data["bus"]
        bus["area"] = area_id_lookup[bus["area"]]
    end
end

"helper function to get all areas IDs"
get_areas_id(data::Dict{String, <:Any})::Vector{Int64} = unique([bus["area"] for (i, bus) in data["bus"]])

"helper function to get all areas IDs"
get_areas_id(data::Dict{Int, <:Any})::Vector{Int64} = collect(keys(data))

"helper function to get all areas IDs"
get_areas_id(pm::AbstractPowerModel)::Vector{Int64} = get_areas_id(pm.data)

"helper function to get the area ID"
get_area_id(data::Dict{String, <:Any})::Int64 = get(data,"area", NaN)

"helper function to get the area ID"
get_area_id(pm::AbstractPowerModel)::Int64 = get_area_id(pm.data)

"helper functions to get the area's local buses"
get_local_bus(data::Dict{String, <:Any}, area::Int64)::Vector{Int64} = [bus["bus_i"] for (i,bus) in data["bus"] if bus["area"] == area]

"helper functions to get the area's local buses"
get_local_bus(pm::AbstractPowerModel, area::Int64)::Vector{Int64} = get_local_bus(pm.data, area)

"helper functions to get the area's neighbor buses"
function get_neighbor_bus(data::Dict{String, <:Any}, local_bus::Vector)::Vector{Int64}
    neighbor_bus = Vector{Int64}()
    for (i,branch) in data["branch"]
        if branch["br_status"] == 1
            if branch["f_bus"] in local_bus && !(branch["t_bus"] in local_bus)
                push!(neighbor_bus,branch["t_bus"])
            elseif !(branch["f_bus"] in local_bus) && branch["t_bus"] in local_bus
                push!(neighbor_bus,branch["f_bus"])
            end
        end
    end
    return neighbor_bus
end

"helper functions to get the area's neighbor buses"
get_neighbor_bus(pm::AbstractPowerModel, local_bus::Vector)::Vector{Int64} = get_neighbor_bus(pm.data, local_bus)

"helper functions to get the area's neighbor buses"
get_neighbor_bus(data::Dict{String, <:Any}, area::Int64)::Vector{Int64} = get_neighbor_bus(data, get_local_bus(data,area))

"helper functions to get the area's neighbor buses"
get_neighbor_bus(pm::AbstractPowerModel, area::Int64)::Vector{Int64} = get_neighbor_bus(pm.data, area)

"helper functions to all areas buses in a dicrionary"
function get_areas_bus(data::Dict{String, <:Any})
    areas_id = get_areas_id(data)
    areas_bus = Dict{Int64, Vector{Int64}}()
    for area in areas_id
        areas_bus[area] = [bus["bus_i"] for (idx,bus) in data["bus"] if bus["area"]==area]
    end
    areas_bus[0] = [bus["bus_i"] for (idx,bus) in data["bus"]]
    return areas_bus
end

"helper functions to all areas buses in a dicrionary"
get_areas_bus(pm::AbstractPowerModel) = get_areas_bus(pm.data)

"get the shared buses and branches between an area and all other areas"
function get_shared_component(data::Dict{String, <:Any}, area_id::Int64)
    areas_id = get_areas_id(data)
    areas_bus = get_areas_bus(data)
    shared_branch = Dict{Int64, Any}()
    shared_bus = Dict{Int64, Any}()
    for area in areas_id
        if area != area_id
            shared_branch[area] = Vector{Int64}(unique([parse(Int64,idx) for (idx,branch) in data["branch"] if branch["br_status"] == 1 && ((branch["f_bus"] in areas_bus[area] && branch["t_bus"] in areas_bus[area_id]) || (branch["f_bus"] in areas_bus[area_id] && branch["t_bus"] in areas_bus[area])) ]))
        else
            shared_branch[area] = Vector{Int64}(unique([parse(Int64,idx) for (idx,branch) in data["branch"] if branch["br_status"] == 1 && xor(branch["f_bus"] in areas_bus[area], branch["t_bus"] in areas_bus[area]) ]))
        end
        
        shared_bus[area] = Vector{Int64}(unique(vcat([branch["f_bus"] for (idx,branch) in data["branch"] if parse(Int64,idx) in shared_branch[area]], [branch["t_bus"] for (idx,branch) in data["branch"] if parse(Int64,idx) in shared_branch[area]] )))
    end
    shared_bus[0] = Vector{Int64}(unique([idx for area in areas_id for idx in shared_bus[area]]))
    shared_branch[0] = Vector{Int64}(unique([idx for area in areas_id for idx in shared_branch[area]]))
    return shared_bus, shared_branch
end

"get the shared buses and branches between defined area and all other areas"
get_shared_component(pm::AbstractPowerModel, area::Int64) = get_shared_component(pm.data, area)

"get the shared buses and branches between defined area and all other areas"
function get_shared_component(data::Dict{String, <:Any})
    area = get_area_id(data)
    get_shared_component(data, area)
end

"get the shared buses and branches between defined area and all other areas"
get_shared_component(pm::AbstractPowerModel) = get_shared_component(pm.data)