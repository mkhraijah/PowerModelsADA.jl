# Data Structure

```@meta
CurrentModule = PowerModelsADA
```


## Input Data 


### Case
The case data is stored in dictionary similar to the one in `PowerModels` with one essential difference. The buses in the data dictionary must contain area key with more than one distinct areas' ID. To load a data file, we use the same `parse_file` function as follow: 

```julia
case_path = "test/data/case14.m"
data = parse_file(case_path)
```

### Partitioning
To check the areas ID in a data dictionary, use `get_areas_id(data)` to get all areas IDs in `data`. If the data dictionary dosen't contain more than one area, we parition the system manually using `assign_area!` or using partitioning algorithm in [KaHyPar.jl](https://github.com/kahypar/KaHyPar.jl) using `partition_system!`. The 

```@docs
assign_area!
partition_system!
```
An example of partition file is shown in [parition example](https://github.com/mkhraijah/PowerModelsADA.jl/blob/main/test/data/case14_2areas.csv).

Before running the distirbuted algorithm, `PowerModelsADA.jl` internally decmpose the original system into subsystems. It decuple the tie-lines by introducing dummy buses and virtual generators at the tie-lines ends. This process is pefromed using `decompose_system` function. 

```@docs
decompose_system
```

## Output Data 
The output of the distributed algorithms is stored in a dictionary where keys are the areas ID and the value are the area data dictionary that contain the results. 


## Iterations Data
To save a specific data during the distributed algorithm (e.g., store the `"shared_variale"` dicrionaryat each iteration), use the option `sace_data::Vector{String}=[]` in the solve function and add the key of the data (e.g., `sace_data=["shared variales"]` ). The output of solve function will contain a dictionary with a key called `"previous_solution"` that constains vectors of the selected saved data ordered by the iteration number.

## Generation Cost
To calculate the objective function of the central algorithm use `calc_dist_gen_cost`.


```@docs
calc_dist_gen_cost
```
To compar the distributed algorithm objective function value with the central OPF, use `compare_solution` to get the absolute relative error. 

```@docs
compare_solution
```