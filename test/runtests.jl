using DistributedPowerModels

import JuMP
import HiGHS
import Ipopt
import PowerModels

using Test

## default setup for solvers
milp_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

const DPM = DistributedPowerModels
PowerModels.silence()

## load test systems
data_14 =  DPM.parse_file("../test/data/case14.m")
data_RTS =  DPM.parse_file("../test/data/case_RTS.m")
data_300 =  DPM.parse_file("../test/data/case300.m")
data_300_3 =  DPM.parse_file("../test/data/case300_3areas.m")

@testset "DistributedPowerModels" begin

## assign area test
    @testset "assign area to busses" begin
        DPM.assign_area!(data_300, "../test/data/case300_8areas.csv")
        area_count = [65, 35, 33, 38, 43, 44, 22, 20]
        test_count = [count(c -> c["area"]  == k, [bus for (i,bus) in data_300["bus"]]) for k in 1:8]
        @test test_count == area_count
    end

## decompose system test
    @testset "decmpose system into subsystems" begin
        @testset "case300" begin
            data_area = DPM.decompose_system(data_300_3)
            bus_size = [169, 83, 69]
            gen_size = [41, 27, 22]
            branch_size = [216, 117, 90]
            load_size = [109, 47, 45]
            neighbor_size = [12, 5, 7]
            test_bus = [length(data_area[i]["bus"]) for i in 1:3]
            test_gen = [length(data_area[i]["gen"]) for i in 1:3]
            test_branch = [length(data_area[i]["branch"]) for i in 1:3]
            test_load = [length(data_area[i]["load"]) for i in 1:3]
            test_neighbor = [length(data_area[i]["neighbor_bus"]) for i in 1:3]

            @test test_bus == test_bus
            @test gen_size == test_gen
            @test branch_size == test_branch
            @test load_size == test_load
            @test neighbor_size == test_neighbor

        end
    end

## paritiotioning test
    @testset "partition system using spectral clustering" begin
        @testset "case_RTS" begin
            DPM.partition_system!(data_RTS, 3, init=[101, 201, 301])

            test_count = [count(c -> c["area"]  == k, [bus for (i,bus) in data_RTS["bus"]]) for k in 1:3]
            test_count_24 = count(==(24),test_count)
            test_count_25 = count(==(25),test_count)
            @test test_count_24 == 2
            @test test_count_25 == 1

        end
    end

## ADMM test
    @testset "admm algorithm with DC power flow" begin

        data_area = DPM.run_dopf_admm(data_14, "DC", milp_solver, 1000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test abs(dist_error) <= 1e0

    end

    @testset "admm algorithm with AC power flow" begin

        data_area = DPM.run_dopf_admm(data_14, "AC", nlp_solver, 1000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "AC", nlp_solver)
        @test abs(dist_error) <= 1e0

    end

    ## ATC test
    @testset "atc algorithm with DC power flow" begin

        data_area = DPM.run_dopf_atc(data_14, "DC", milp_solver, 1.1, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test abs(dist_error) <= 1e0

    end

    @testset "atc algorithm with SOC relaxation of power flow" begin

        data_area = DPM.run_dopf_atc(data_14, "SOCP", nlp_solver, 1.1, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "SOCP", nlp_solver)
        @test abs(dist_error) <= 1e-1

    end

    ## APP test
    @testset "app algorithm with DC power flow" begin

        data_area = DPM.run_dopf_app(data_14, "DC", milp_solver, 1000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test abs(dist_error) <= 1e0

    end

    @testset "app algorithm with DC power flow" begin

        data_area = DPM.run_dopf_app(data_14, "ACR", nlp_solver, 1000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "ACR", nlp_solver)
        @test abs(dist_error) <= 1e-1

    end
end
