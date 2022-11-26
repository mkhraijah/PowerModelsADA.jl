###############################################################################
#                   Methods for sharing data between areas                    #
###############################################################################

"prepare the shared data with or without serialization"
function prepare_shared_data(data::Dict{String, <:Any}, to_area::Int64; serialize::Bool=false)
    shared_data_key = filter(x -> startswith(string(x), "shared"), keys(data))
    shared_data = Dict{String, Any}([key => get(data[key], string(to_area), Dict()) for key in shared_data_key])

    if serialize
        shared_data = _serialize_shared_data!(shared_data)
    end

    return shared_data
end

function _serialize_shared_data!(shared_data::Dict{String, <:Any})
    io = IOBuffer()
    Serialization.serialize(io, shared_data)
    return take!(io)
end

"deserialize and store the received data in the local data"
function receive_shared_data!(data::Dict{String, <:Any}, shared_data::Vector, from_area::Int64)
    shared_data = Serialization.deserialize(IOBuffer(shared_data))
    receive_shared_data!(from_area, shared_data, data)
end

"store received data in the local data dictionary"
function receive_shared_data!(data::Dict{String, <:Any}, shared_data::Dict{String, <:Any}, from_area::Int64)
    for shared_data_key in keys(shared_data)
        key = "received_$(shared_data_key[8:end])"
        if haskey(data, key)
            data[key][string(from_area)] = shared_data[shared_data_key]
        end
    end
end