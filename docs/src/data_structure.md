# Data Structure

```@meta
CurrentModule = PowerModelsADA
```


## Input Data 


### Case
`PowerModelsADA` uses a dictionary of dictionaries to store the case data and the subproblem data. The case data dictionary is similar to the one in `PowerModels` with area assignment for each bus. The buses in the data dictionary must contain an area key with more than one distinct area ID. The subproblem data is similar to the case data with additional information. The area data contains the area-specific data and the distributed algorithm parameters. 

To load a data file, we use `parse_file` function as follow: 

```julia
case_path = "test/data/case14.m"
data = parse_file(case_path)
```

### Partitioning
To check the areas ID in a data dictionary, use `get_areas_id(data)` to get all areas' IDs in `data`. If the data dictionary doesn't contain more than one area, we can partition the system manually using `assign_area!` function or using the partitioning algorithm in [KaHyPar.jl](https://github.com/kahypar/KaHyPar.jl) via `partition_system!` function. An example of partition file is shown in [partition example](https://github.com/mkhraijah/PowerModelsADA.jl/blob/main/test/data/case14_2areas.csv).

```@docs
assign_area!
partition_system!
```

Before running the distributed algorithm, `PowerModelsADA` internally decomposes the original system into subsystems using `decompose_system` function. the function decouples the tie-lines between two areas by introducing dummy buses and virtual generators at the tie-lines' ends.

```@docs
decompose_system
```

## Output Data 
The output of the distributed algorithms is stored in a dictionary. The dictionary's keys are the areas ID, and the dictionary's values are the areas data dictionary with the results stored in `solution` dictionary. 


## Saving Iterations Data
To save a specific data during the distributed algorithm (e.g., store the `"shared_variable"` dictionary each iteration), use the option `save_data::Vector{String}=[]` in the solve function and add the key of the data (e.g., `save_data=["shared variable"]`). The output of the solve function will contain a dictionary with a key called `"previous_solution"` that contains vectors of the selected saved data ordered by the iteration number.

## Generation Cost
To calculate the objective function of the central algorithm use `calc_dist_gen_cost`.


```@docs
calc_dist_gen_cost
```
To compare the distributed algorithm objective function value with the central OPF, use `compare_solution` to get the absolute value of the relative error. 

```@docs
compare_solution
```

