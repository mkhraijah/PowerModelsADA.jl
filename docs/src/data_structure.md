# Data Structure

```@meta
CurrentModule = PMADA
```


## Input Data 


### Case
The input data is the same as `PowerModels` with one essential exeption. The buses in the `PowerModels` dictionary must contain area key with more than one area number. To load a data file, we use the same `PowerModels` methods. 
The `PowerModels` data can be loaded as follow: 

```julia
case_path = "test/data/case14.m"
data = parse_file(case_path)
```

### Partitioning
To check the areas in a `PowerModels` data, use `PMADA.get_areas_id(data)` to get all areas ids in `data`. If the data dictionary dosen't contain more than one area, there are two methods to parition the system `assign_area!` or `partition_system!`

```@docs
assign_area!
partition_system!
```
An example of partition file is shown in [parition example](https://github.com/mkhraijah/PMADA.jl/blob/main/test/data/case14_2areas.csv).

Before running the distirbuted algorithm, `PMADA.jl` internally decmpose the original system into subsystems. It decuple the tie-lines by introducing dummy buses and virtual generators at the tie-lines ends. This process is pefromed using `decompose_system` function. 

```@docs
decompose_system
```

## Output Data 
The output of the distributed algorithms is stored in a dictionary the keys are the area id and the value are the area data dictionary that contain the results. 


## Historical Data
TODO

## Generation Cost
To calculate the objective function of the central algorithm use `calc_dist_gen_cost`.


```@docs
calc_dist_gen_cost
```
To compar the distributed algorithm objective function value with the central OPF, use `compare_solution` to get the absolute relative error. 

```@docs
compare_solution
```