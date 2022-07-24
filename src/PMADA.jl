
module PMADA


import JuMP
import PowerModels as _PM
import Serialization
import LinearAlgebra
import DelimitedFiles
import Clustering: kmeans

import PowerModels: AbstractPowerModel, DCPPowerModel, ACPPowerModel, ACRPowerModel, SOCWRPowerModel, QCRMPowerModel, SDPWRMPowerModel

import PowerModels: parse_file, ids


include("core/base.jl")
include("core/variables.jl")
include("core/objective.jl")
include("core/data.jl")
include("core/data_sharing.jl")
include("core/util.jl")

include("algorithms/admm_methods.jl")
include("algorithms/atc_methods.jl")
include("algorithms/app_methods.jl")


"Suppresses information and warning messages output by PowerModels"
function silence!()
    @warn "Suppressing information and warning messages output by PowerModels and PowerModelsDistribution for the rest of this session."
    _PM.silence()
end


end
