# PowerModelsAlternatingDistributedAlgorithms.jl (PMADA.jl)

A library to run different distributed algorithms to solve optimal power flow using PowerModels and JuMP. Currently, the library uses three distributed algorithms: Alternating Direction Method of Multipliers (ADMM), Auxiliary Problem Principle (APP), and Analytical Target Cascading (ATC). 


Status:
[![CI](https://github.com/mkhraijah/DistributedPowerModels/workflows/CI/badge.svg)](https://github.com/mkhraijah/DistributedPowerModels/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/mkhraijah/DistributedPowerModels/branch/main/graph/badge.svg)](https://codecov.io/gh/mkhraijah/DistributedPowerModels)
</p>

## Dependencies
* InfrastructureModels v0.7.4
* PowerModels v0.19.5
* JuMP v1.0.0
* LinearAlgebra

## Usage

There are two modes of running the distributed algorithms: global and local modes. The global mode runs the distributed algorithms by passing the information of all subsystems, while the local mode assumes each subsystem is handled by a separate object and can be distributed over different threads or processors. 

The distirbuted algorithms can be run in global mode as follow: 

`run_dopf_admm(data, model_type, optimizer, alpha=1000, tol=1e-4, max_iteration=1000, verbose = true)`


## Configuration

To configure the model, change the directory in Julia to the library folder and activate the environment using `Pkg. activate(“.”)`, and load the module using `include("PMADA.jl")`. The case study needs to be parsed using `parse_file(casename)`. 

*Note: the parsed case study needs to be divided into subareas using the key “area” of each bus. You can use the helper function `assign_area!(data, partition_file)` to load a partition csv file with bus and area pairs in each row and assign the bus area accordingly.  

### Power Flow formulation 

The power flow formulation is selected using the variable ` pf_model`. Currently, the code supports the following models from PowerModels: 
`ACPPowerModel`, `ACRPowerModel`, `DCPPowerModel`, `SOCWRPowerModel`, `QCRMPowerModel`, `SDPWRMPowerModel`. See https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/core/types.jl for more details about the supported power flow formulations. 


## Examples 

For examples on how to run the code see `example/example.jl`.
