## Methods for data wrangling, sharing and obtaining an information from data


## wrapper for PowerModel.parse_file
parse_file(case_path::String) = _PM.parse_file(case_path)

## Assign area to the PowerModel data using a dictionary with (bus => area) Int pairs
function assign_area!(data::Dict, partition::Dict)
    for i in keys(data["bus"])
        data["bus"][i]["area"] = partition[parse(Int64,i)]
    end
end

## Assign area to the PowerModel data using a CVS file with buses and area id
function assign_area!(data::Dict, partition_path::String)
    partition = DelimitedFiles.readdlm(partition_path, ',', Int, '\n')
    assign_area!(data, partition)
end

## Assign area to the PowerModel data using a vector with (bus => area) pairs
function assign_area!(data::Dict, partition::Vector{Pair{Int64, Int64}})
    assign_area!(data, Dict(partition))
end

## Assign area to the PowerModel data using a matrix with [bus, area] colmuns or rows
function assign_area!(data::Dict, partition::Array{Int64, 2})
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
function decompose_system(data::Dict{String, <:Any}, area::Int)

    # idintify local buses
    local_bus = get_local_bus(data, area)
    neighbor_bus = get_neighbor_bus(data, area)

    ## generators list
    virtual_gen =  add_virtual_gen!(data,neighbor_bus, area)
    area_gen = Dict{String, Any}([i => gen for (i,gen) in data["gen"] if gen["gen_bus"] in local_bus])

    data_area = Dict{String,Any}()
    data_area["area"] = area
    data_area["name"]= "$(data["name"])_area_$area"
    data_area["source_version"] = data["source_version"]
    data_area["source_type"] = data["source_type"]
    data_area["baseMVA"] = data["baseMVA"]
    data_area["per_unit"] = data["per_unit"]
    data_area["bus"] = Dict{String,Any}([j => bus for (j,bus) in data["bus"] if bus["bus_i"] in [local_bus;neighbor_bus]])
    data_area["branch"] = Dict{String,Any}([j => branch for (j,branch) in data["branch"] if branch["f_bus"] in local_bus || branch["t_bus"] in local_bus])
    data_area["gen"] = merge(area_gen,virtual_gen)
    data_area["shunt"] = Dict{String, Any}([i => shunt for (i,shunt) in data["shunt"] if shunt["shunt_bus"] in local_bus])
    data_area["load"] = Dict{String, Any}([i => load for (i,load) in data["load"] if load["load_bus"] in local_bus])
    data_area["storage"]= Dict{String, Any}([i => storage for (i,storage) in data["storage"] if gen["storage_bus"] in local_bus])
    data_area["switch"]=Dict{String, Any}([i => switch for (i,switch) in data["switch"] if gen["switch_bus"] in local_bus])
    data_area["dcline"]= Dict{String, Any}([i => dcline for (i,dcline) in data["dcline"] if dcline["f_bus"] in local_bus || dcline["t_bus"] in local_bus ] )

    data_area["neighbor_bus"] = Dict{String,Any}([j => bus for (j,bus) in data["bus"] if bus["bus_i"] in neighbor_bus])

    return data_area
end

# add virtual geneartors at the neighboring buses of an area
function add_virtual_gen!(data::Dict{String, <:Any},neighbor_bus::Vector, area::Int)
    max_gen_ind = maximum([parse(Int,i) for i in keys(data["gen"])])
    virtual_gen = Dict{String, Any}()
    cost_model = data["gen"]["$max_gen_ind"]["model"]
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

## methot to get the shared data with or without serialization
function send_shared_data(from::Int64, to::Int64, pm::AbstractPowerModel; serialize::Bool=false)
    shared_data = Dict{Symbol, Any}();
    for j in keys(pm.ext[:primal_shared_variable][to]) # loop through variables (va, vm, p, q)
        shared_data[j] = Dict{Any, Any}();
        for k in keys(pm.ext[:primal_shared_variable][to][j]) # loop through the shared variables with area "to"
                shared_data[j][k] = pm.ext[:primal_shared_variable][from][j][k]
        end
    end
    if serialize
        # IObuffer function to convert object to byte streams
        io = IOBuffer()
        # Serialize function takes stream and value as parameters
        serialize(io, shared_data)
        # take! Function fetches IOBUffer contents as Byte array
        shared_data = take!(io)
    end
    return shared_data
end

## Method to deserialize and store the received data in the local pm object
function receive_shared_data!(from::Int64, shared_data::Vector, pm::AbstractPowerModel)
    shared_data = deserialize(IOBuffer(shared_data))
    receive_shared_data!(from,shared_data, pm)
end

function receive_shared_data!(from::Int64, shared_data::Dict, pm::AbstractPowerModel)
    for i in keys(pm.ext[:primal_shared_variable][from])
        for j in keys(pm.ext[:primal_shared_variable][from][i])
            if !isnan(shared_data[i][j])
                pm.ext[:primal_shared_variable][from][i][j] = shared_data[i][j]
            end
        end
    end
end


## helper functions to handle area ids, local buses, neighbor buses
get_areas_id(data::Dict) = unique([bus["area"] for (i, bus) in data["bus"]])

get_areas_id(pm::AbstractPowerModel) = unique([bus["area"] for (i, bus) in pm.data["bus"]])

get_area_id(data::Dict) = get(data,"area", NaN)

get_area_id(pm::AbstractPowerModel) = get(pm.data,"area", NaN)

get_local_bus(data::Dict,area::Int) = [bus["bus_i"] for (i,bus) in data["bus"] if bus["area"] == area]

get_local_bus(pm::AbstractPowerModel,area::Int) = [bus["bus_i"] for (i,bus) in pm.data["bus"] if bus["area"] == area]

function get_neighbor_bus(data::Dict, local_bus::Vector)
    neighbor_bus = []
    for (i,branch) in data["branch"]
        if branch["f_bus"] in local_bus && !(branch["t_bus"] in local_bus)
            push!(neighbor_bus,branch["t_bus"])
        elseif !(branch["f_bus"] in local_bus) && branch["t_bus"] in local_bus
            push!(neighbor_bus,branch["f_bus"])
        end
    end
    return neighbor_bus
end

get_neighbor_bus(pm::AbstractPowerModel, area::Int) = get_neighbor_bus(pm.data, area)

get_neighbor_bus(data::Dict, area::Int) = get_neighbor_bus(data, get_local_bus(data,area))

function get_areas_bus(data::Dict)
    areas_id = get_areas_id(data)
    areas_bus = Dict{Int64, Any}()
    for i in areas_id
        areas_bus[i] = [bus["bus_i"] for (j,bus) in data["bus"] if bus["area"]==i]
    end
    return areas_bus
end

get_areas_bus(pm::AbstractPowerModel) = get_areas_bus(pm.data)

## Method to get the shared buses and branches between defined area and all other areas in pm.data
function get_shared_component(pm::AbstractPowerModel, area::Int64)
    areas_id = get_areas_id(pm)
    areas_bus = get_areas_bus(pm)
    shared_branch = Dict{Int64, Any}()
    shared_bus = Dict{Int64, Any}()
    for i in areas_id
        if i != area
            shared_branch[i] = unique([parse(Int64,j) for (j,branch) in pm.data["branch"] if (branch["f_bus"] in areas_bus[i] && branch["t_bus"] in areas_bus[area]) || (branch["f_bus"] in areas_bus[area] && branch["t_bus"] in areas_bus[i]) ])
        else
            shared_branch[i] = unique([parse(Int64,j) for (j,branch) in pm.data["branch"] if xor(branch["f_bus"] in areas_bus[i], branch["t_bus"] in areas_bus[i]) ])
        end
            shared_bus[i] = unique(vcat([branch["f_bus"] for (j,branch) in pm.data["branch"] if parse(Int64,j) in shared_branch[i]], [branch["t_bus"] for (j,branch) in pm.data["branch"] if parse(Int64,j) in shared_branch[i]] ))
    end
    return shared_bus, shared_branch
end

function get_shared_component(pm::AbstractPowerModel)
    area = get_area_id(pm)
    get_shared_component(pm, area)
end

## helper method to find the correct power flow model
function pf_formulation(pf::DataType)
    if pf <: AbstractPowerModel
        pf
    else
        error("Power Flow model is not identified")
    end
end

function pf_formulation(pf::String)
    if pf == "DC"
        _PM.DCPPowerModel
    elseif pf == "AC" || pf == "ACP"
        _PM.ACPPowerModel
    elseif pf == "ACR"
        _PM.ACRPowerModel
    elseif pf == "SOC"
        _PM.SOCWRPowerModel
    elseif pf == "QC"
        _PM.QCRMPowerModel
    elseif pf == "SDP"
        _PM.SDPWRMPowerModel
    else
        error("Power Flow model is not identified")
    end
end
