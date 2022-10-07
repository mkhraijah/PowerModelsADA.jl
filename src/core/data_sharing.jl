###############################################################################
#                   Methods for sharing data between areas                    #
###############################################################################

"get the shared data with or without serialization"
function send_shared_data(from::Int64, to::Int64, data::Dict{String, <:Any}; serialize::Bool=false)

    shared_data_key = filter(x -> startswith(string(x), "shared"), keys(data))

    shared_data = Dict{String, Any}([key => Dict{String, Any}([
        variable => Dict{String, Any}([
            idx => data[key][string(from)][variable][idx] 
            for idx in keys(data[key][string(to)][variable]) ])
        for variable in keys(data[key][string(to)]) ])
    for key in shared_data_key ])

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

"deserialize and store the received data in the local data"
function receive_shared_data!(from::Int64, shared_data::Vector, data::Dict{String, <:Any})
    shared_data = Serialization.deserialize(IOBuffer(shared_data))
    receive_shared_data!(from,shared_data, data)
end

"store received data in the local data dictionary"
function receive_shared_data!(from::Int64, shared_data::Dict{String, <:Any}, data::Dict{String, <:Any})
    for shared_data_key in keys(shared_data)
        for variable in keys(data[shared_data_key][string(from)])
            for idx in keys(data[shared_data_key][string(from)][variable])
                if !isnan(shared_data[shared_data_key][variable][idx])
                    if isa(shared_data[shared_data_key][variable][idx], String)
                        value = parse(Float64, shared_data[shared_data_key][variable][idx])
                    else
                        value = shared_data[shared_data_key][variable][idx]
                    end
                    data[shared_data_key][string(from)][variable][idx] = value
                end
            end
        end
    end
end