###############################################################################
#                   Methods for sharing data between areas                    #
###############################################################################

## methot to get the shared data with or without serialization
function send_shared_data(from::Int64, to::Int64, data::Dict{String, <:Any}; serialize::Bool=false)

    shared_data = Dict{String, Any}()
    for j in keys(data["shared_primal"][string(to)]) # loop through variables (va, vm, p, q)
        shared_data[j] = Dict{String, Any}()
        for k in keys(data["shared_primal"][string(to)][j]) # loop through the shared variables with area "to"
                shared_data[j][k] = data["shared_primal"][string(from)][j][k]
        end
    end

    if serialize
        # IObuffer function to convert object to byte streams
        io = IOBuffer()
        # Serialize function takes stream and value as parameters
        Serialization.serialize(io, shared_data)
        # take! Function fetches IOBUffer contents as Byte array
        shared_data = take!(io)
    end

    return shared_data
end

## Method to deserialize and store the received data in the local pm object
function receive_shared_data!(from::Int64, shared_data::Vector, data::Dict{String, <:Any})
    shared_data = Serialization.deserialize(IOBuffer(shared_data))
    receive_shared_data!(from,shared_data, data)
end

function receive_shared_data!(from::Int64, shared_data::Dict, data::Dict{String, <:Any})
    for comp in keys(data["shared_primal"][string(from)])
        for ids in keys(data["shared_primal"][string(from)][comp])
            for vstring in keys(data["shared_primal"][string(from)][comp][ids])
                if !isnan(shared_data[comp][ids][vstring])
                    data["shared_primal"][string(from)][comp][ids][vstring] = shared_data[comp][ids][vstring]
                end
            end
        end
    end
end

## Method to update the primal variable after solving the subproblem
function update_primal_shared!(data::Dict)
    area_id = string(get_area_id(data))
    for comp in keys(data["shared_primal"][area_id])
        for ids in keys(data["shared_primal"][area_id][comp])
            for vstring in keys(data["shared_primal"][area_id][comp][ids])
                data["shared_primal"][area_id][comp][ids][vstring] = data[comp][ids][vstring]
            end
        end
    end
end
