###############################################################################
# Methods for data wrangling, sharing and obtaining an information from data  #
###############################################################################


## assign area to the PowerModel data using a dictionary with (bus => area) Int pairs
function assign_area!(data::Dict{String, <:Any}, partition::Dict)
    for i in keys(data["bus"])
        data["bus"][i]["area"] = partition[parse(Int64,i)]
    end
end

## assign area to the PowerModel data using a CVS file with buses and area id
function assign_area!(data::Dict{String, <:Any}, partition_path::String)
    partition = DelimitedFiles.readdlm(partition_path, ',', Int, '\n')
    assign_area!(data, partition)
end


## Assign area to the PowerModel data using a vector with (bus => area) pairs
function assign_area!(data::Dict{String, <:Any}, partition::Vector{Pair{Int64, Int64}})
    assign_area!(data, Dict(partition))
end

## Assign area to the PowerModel data using a matrix with [bus, area] colmuns or rows
function assign_area!(data::Dict{String, <:Any}, partition::Array{Int64, 2})
    if size(partition)[2] != 2 && length(data["bus"]) != 2
        partition = partition'
        if size(partition)[2] != 2
            #through error
            error("Partitioning data doesn't contin correct area assignment")
        end
    end
    assign_area!(data, Dict(partition[i,1] => partition[i,2] for i in 1:size(partition)[1] ))
end

## method to decompose a subsystem with area id
function decompose_system(data::Dict{String, <:Any}, area_id::Int)

    # idintify local buses
    local_bus = get_local_bus(data, area_id)
    neighbor_bus = get_neighbor_bus(data, area_id)

    ## generators list
    virtual_gen =  add_virtual_gen!(data,neighbor_bus, area_id)
    area_gen = Dict{String, Any}([i => gen for (i,gen) in data["gen"] if gen["gen_bus"] in local_bus])

    data_area = Dict{String,Any}()
    data_area["area"] = area_id
    data_area["name"]= "$(data["name"])_area_$area_id"
    data_area["source_version"] = data["source_version"]
    data_area["source_type"] = data["source_type"]
    data_area["baseMVA"] = data["baseMVA"]
    data_area["per_unit"] = data["per_unit"]
    data_area["bus"] = Dict{String,Any}([j => bus for (j,bus) in data["bus"] if bus["bus_i"] in [local_bus;neighbor_bus]])
    data_area["branch"] = Dict{String,Any}([j => branch for (j,branch) in data["branch"] if branch["f_bus"] in local_bus || branch["t_bus"] in local_bus])
    data_area["gen"] = merge(area_gen, virtual_gen)
    data_area["shunt"] = Dict{String, Any}([i => shunt for (i,shunt) in data["shunt"] if shunt["shunt_bus"] in local_bus])
    data_area["load"] = Dict{String, Any}([i => load for (i,load) in data["load"] if load["load_bus"] in local_bus])
    data_area["storage"]= Dict{String, Any}([i => storage for (i,storage) in data["storage"] if gen["storage_bus"] in local_bus])
    data_area["switch"]=Dict{String, Any}([i => switch for (i,switch) in data["switch"] if gen["switch_bus"] in local_bus])
    data_area["dcline"]= Dict{String, Any}([i => dcline for (i,dcline) in data["dcline"] if dcline["f_bus"] in local_bus || dcline["t_bus"] in local_bus ] )

    data_area["local_bus"] = local_bus
    data_area["neighbor_bus"] = neighbor_bus

    return data_area
end

# add virtual geneartors at the neighboring buses of an area
function add_virtual_gen!(data::Dict{String, <:Any},neighbor_bus::Vector, area_id::Int)
    max_gen_ind = maximum([parse(Int,i) for i in keys(data["gen"])])
    virtual_gen = Dict{String, Any}()
    cost_model = data["gen"][string(max_gen_ind)]["model"]
    max_flow = 10*sum(load["pd"] for (i,load) in data["load"])
    if cost_model == 1
        for i in neighbor_bus
            virtual_gen[string(i+max_gen_ind)] = Dict{String, Any}("ncost" => 2, "qc1max" => 0.0, "pg" => 0, "model" => cost_model, "shutdown" => 0.0, "startup" => 0.0, "qc2max" => 0.0, "ramp_agc" => 0.0, "qg" => 0.0, "gen_bus" => i, "pmax" => max_flow, "ramp_10" => 0.0, "vg" => 1.05, "mbase" => data["baseMVA"], "source_id" => Any["gen", max_gen_ind+i], "pc2" => 0.0, "index" => i+max_gen_ind, "cost" => [0.0; 0.0; 0.0; 0.0], "qmax" => max_flow, "gen_status" => 1, "qmin" => -max_flow, "qc1min" => 0.0, "qc2min" => 0.0, "pc1" => 0.0, "ramp_q" => 0.0, "ramp_30" => 0.0, "pmin" => -max_flow, "apf" => 0.0)
        end
    else
        for i in neighbor_bus
            virtual_gen[string(i+max_gen_ind)] = Dict{String, Any}("ncost" => 3, "qc1max" => 0.0, "pg" => 0, "model" => cost_model, "shutdown" => 0.0, "startup" => 0.0, "qc2max" => 0.0, "ramp_agc" => 0.0, "qg" => 0.0, "gen_bus" => i, "pmax" => max_flow, "ramp_10" => 0.0, "vg" => 1.05, "mbase" => data["baseMVA"], "source_id" => Any["gen", max_gen_ind+i], "pc2" => 0.0, "index" => i+max_gen_ind, "cost" => [0.0; 0.0; 0.0], "qmax" => max_flow, "gen_status" => 1, "qmin" => -max_flow, "qc1min" => 0.0, "qc2min" => 0.0, "pc1" => 0.0, "ramp_q" => 0.0, "ramp_30" => 0.0, "pmin" => -max_flow, "apf" => 0.0)
        end
    end
    return virtual_gen
end
