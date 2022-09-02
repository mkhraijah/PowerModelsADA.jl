using PMADA

import JuMP
import HiGHS
import Ipopt
import PowerModels

using Test

## default setup for solvers
milp_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

PowerModels.silence()

## load test systems
data_14 =  PMADA.parse_file("../test/data/case14.m")
data_RTS =  PMADA.parse_file("../test/data/case_RTS.m")

@testset "PMADA" begin

## assign area test
    @testset "assign area to busses" begin
        PMADA.assign_area!(data_14, "../test/data/case14_2areas.csv")
        area_bus = PMADA.get_areas_bus(data_14)
        area_1 = [1, 2, 3, 4, 5]
        area_2 = [6, 7, 8, 9, 10, 11, 12, 13, 14]
        @test sort(area_bus[1]) == area_1
        @test sort(area_bus[2]) == area_2
    end

## decompose system test 
    @testset "decmpose system into subsystems" begin
        @testset "data_RTS" begin
            data_area = PMADA.decompose_system(data_RTS)
            bus_size = [28, 28, 27]
            gen_size = [56, 41, 71]
            branch_size = [42, 42, 41]
            load_size = [17, 17, 17]
            neighbor_size = [4, 4, 2]
            test_bus = [length(data_area[i]["bus"]) for i in 1:3]
            test_gen = [length(data_area[i]["gen"]) for i in 1:3]
            test_branch = [length(data_area[i]["branch"]) for i in 1:3]
            test_load = [length(data_area[i]["load"]) for i in 1:3]
            @test test_bus == test_bus
            @test gen_size == test_gen
            @test branch_size == test_branch
            @test load_size == test_load
        end
    end

## paritiotioning test
    @testset "partition system" begin
        @testset "case_RTS" begin
            PMADA.partition_system!(data_RTS, 3)
            test_count = [count(c -> c["area"]  == k, [bus for (i,bus) in data_RTS["bus"]]) for k in 1:3]
            test_count_24 = count(==(24),test_count)
            test_count_25 = count(==(25),test_count)
            @test test_count_24 == 2
            @test test_count_25 == 1
        end
    end

## ADMM test
    @testset "admm algorithm with DC power flow" begin
        data_area = PMADA.solve_dopf_admm(data_14, PMADA.DCPPowerModel, milp_solver; alpha=1000, tol=1e-3, max_iteration=1000, verbose = false)
        dist_cost = PMADA.calc_dist_gen_cost(data_area)
        @test  isapprox(dist_cost, 7642.59, atol =5)
    end

    @testset "admm algorithm with AC power flow" begin
        data_coordinator, data_area = PMADA.solve_dopf_admm_coordinated(data_14, PMADA.ACPPowerModel, nlp_solver; alpha=1000, tol=1e-3, max_iteration=1000, verbose = false)
        dist_cost = PMADA.calc_dist_gen_cost(data_area)
        @test isapprox(dist_cost, 8081.52, atol =5)
    end

    ## ATC test
    @testset "atc algorithm with DC power flow" begin
        data_area = PMADA.solve_dopf_atc(data_14, PMADA.DCPPowerModel, milp_solver; alpha=1.1, tol=1e-3, max_iteration=1000, verbose = false)
        dist_cost = PMADA.calc_dist_gen_cost(data_area)
        @test isapprox(dist_cost, 7642.59, atol =5)
    end

    @testset "atc algorithm with SOC relaxation of power flow" begin
        data_coordinator, data_area = PMADA.solve_dopf_atc_coordinated(data_14, PMADA.SOCWRPowerModel, nlp_solver; alpha=1.1, tol=1e-3, max_iteration=1000, verbose = false)
        dist_cost = PMADA.calc_dist_gen_cost(data_area)
        @test isapprox(dist_cost, 8075.12, atol =5)
    end

    ## APP test
    @testset "app algorithm with DC power flow" begin
        data_area = PMADA.solve_dopf_app(data_14, PMADA.DCPPowerModel, milp_solver; alpha=1000, tol=1e-3, max_iteration=1000, verbose = false)
        dist_cost = PMADA.calc_dist_gen_cost(data_area)
        @test isapprox(dist_cost, 7642.59, atol =5)
    end

end
