###############################################################################
#                   Methods for sharing data between areas                    #
###############################################################################

## methot to get the shared data with or without serialization
function send_shared_data(from::Int64, to::Int64, data::Dict{String, <:Any}; serialize::Bool=false)
    shared_data = Dict{Symbol, Any}();
    for j in keys(data["shared_primal"][to]) # loop through variables (va, vm, p, q)
        shared_data[j] = Dict{Any, Any}();
        for k in keys(data["shared_primal"][to][j]) # loop through the shared variables with area "to"
                shared_data[j][k] = data["shared_primal"][from][j][k]
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
    for i in keys(data["shared_primal"][from])
        for j in keys(data["shared_primal"][from][i])
            if !isnan(shared_data[i][j])
                data["shared_primal"][from][i][j] = shared_data[i][j]
            end
        end
    end
end
