using DistributedPowerModels

import HiGHS
import Ipopt


using Test

## default setup for solvers
milp_solver = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag"=>false)
nlp_solver = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

## load test systems
data_WB5 = DPM.parse_file("test/data/WB5.m")
data_14 =  DPM.parse_file("test/data/case14.m")
data_30 =  DPM.parse_file("test/data/case30.m")
data_RTS =  DPM.parse_file("test/data/case_RTS.m")
data_118 =  DPM.parse_file("test/data/case118.m")
data_118_3 = DPM.parse_file("test/data/case118_3areas.m")
data_300 =  DPM.parse_file("test/data/case300.m")
data_300_3 =  DPM.parse_file("test/data/case300_3areas.m")


@testset "DistributedPowerModels" begin

## assign area test
    @testset "assign area to busses for case300" begin
        DPM.assign_area!(data_300, "test/data/case300_8areas.csv")
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

        @testset "case118" begin
            data_area = DPM.decompose_system(data_118_3)
            bus_size = [39, 45, 48]
            gen_size = [19, 23, 26]
            branch_size = [52, 69, 73]
            load_size = [29, 30, 40]
            neighbor_size = [4, 7, 5]
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
    @testset "parition system using spectural clustering" begin
        @testset "case_RTS" begin
            DPM.partition_system!(data_RTS, 3)

            test_count = [count(c -> c["area"]  == k, [bus for (i,bus) in data_RTS["bus"]]) for k in 1:3]
            test_count_24 = count(==(24),test_count)
            test_count_25 = count(==(25),test_count)
            @test test_count_24 == 2
            @test test_count_25 == 1

        end
    end

## sharing data
    @testset "sharing data between subsystems" begin
        data_area = DPM.decompose_system(data_RTS)
        for i in keys(data_area)
            DPM.initialize_dpm!(data_area[i], "DC")
        end
        pm = DPM.update_subproblem(data_area[1], "DC", DPM.build_dopf_admm)
        DPM.solve_subproblem!(pm, milp_solver)
        DPM.update_solution!(data_area[i], pm)

        shared_data = DPM.send_shared_data(1, 2, data_area[1], serialize = false)
        va = [0.0, 0.1644, 0.1753, 0.1675, -0.0481, 0.3257]
        p = [-2.1952, -1.3, -2.1401]

        va_test = collect(values(shared_data[:va]))
        p_test = collect(values(shared_data[:p]))
        @test isapprox(va, va_test, atol=1e-2)
        @test isapprox(p, p_test, atol=1e-2)
        shared_data = DPM.send_shared_data(1, 2, data_area[1], serialize = true)
        DPM.receive_shared_data!(1, shared_data, data_area[2])
        va_test = collect(values(data_area[2]["shared_primal"][1][:va]))
        p_test =  collect(values(data_area[2]["shared_primal"][1][:p]))
        @test isapprox(va, va_test, atol=1e-2)
        @test isapprox(p, p_test, atol=1e-2)

    end


## ADMM test
    @testset "admm algorithm with DC power flow" begin
        data_area = DPM.run_dopf_admm(data_WB5, "DC", milp_solver, 100, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_admm(data_14, "DC", milp_solver, 1000, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_admm(data_30, "DC", milp_solver, 1000, tol = 1e-6, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "DC", milp_solver)
        @test dist_error <= 1e-4

        data_area = DPM.run_dopf_admm(data_118, "DC", milp_solver, 100000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "DC", milp_solver)
        @test dist_error <= 1e0


    end

    @testset "admm algorithm with AC power flow" begin
        data_area = DPM.run_dopf_admm(data_WB5, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_admm(data_14, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_30, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_118, "AC", nlp_solver, 100000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "AC", nlp_solver)
        @test dist_error <= 1e0

    end

    @testset "admm algorithm with SOC relaxation of power flow" begin
        data_area = DPM.run_dopf_admm(data_WB5, "SOCP", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_14, "SOCP", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_30, "SOCP", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-0

        data_area = DPM.run_dopf_admm(data_118, "SOCP", nlp_solver, 100000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e0

    end

    @testset "admm algorithm with QC relaxation of power flow" begin
        data_area = DPM.run_dopf_admm(data_WB5, "QC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_14, "QC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_admm(data_30, "QC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-0
    end


    ## ATC test
    @testset "atc algorithm with DC power flow" begin
        data_area = DPM.run_dopf_atc(data_WB5, "DC", milp_solver, 1.1, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_atc(data_14, "DC", milp_solver, 1.1, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_atc(data_30, "DC", milp_solver, 1.05, tol = 1e-6, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "DC", milp_solver)
        @test dist_error <= 1e-4

        data_area = DPM.run_dopf_atc(data_118, "DC", milp_solver, 1.05, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "DC", milp_solver)
        @test dist_error <= 1e0

    end

    @testset "atc algorithm with AC power flow" begin
        data_area = DPM.run_dopf_atc(data_WB5, "AC", nlp_solver, 1.1, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-4

        data_area = DPM.run_dopf_atc(data_14, "AC", nlp_solver, 1.05, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_30, "AC", nlp_solver, 1.02, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_118, "AC", nlp_solver, 1.02, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "AC", nlp_solver)
        @test dist_error <= 1e0

    end

    @testset "atc algorithm with SOC relaxation of power flow" begin
        data_area = DPM.run_dopf_atc(data_WB5, "SOCP", nlp_solver, 1.1, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_14, "SOCP", nlp_solver, 1.05, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_30, "SOCP", nlp_solver, 1.02, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-1

        data_area = DPM.run_dopf_atc(data_118, "SOCP", nlp_solver, 1.02, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

    end

    @testset "atc algorithm with QC relaxation of power flow" begin
        data_area = DPM.run_dopf_atc(data_WB5, "QC", nlp_solver, 1.02, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_14, "QC", nlp_solver, 1.01, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_atc(data_30, "QC", nlp_solver, 1.02, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2
    end



    ## APP test
    @testset "app algorithm with DC power flow" begin
        data_area = DPM.run_dopf_app(data_WB5, "DC", milp_solver, 100, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_app(data_14, "DC", milp_solver, 1000, tol = 1e-8, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "DC", milp_solver)
        @test dist_error <= 1e-6

        data_area = DPM.run_dopf_app(data_30, "DC", milp_solver, 1000, tol = 1e-6, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "DC", milp_solver)
        @test dist_error <= 1e-4

        data_area = DPM.run_dopf_app(data_118, "DC", milp_solver, 100000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "DC", milp_solver)
        @test dist_error <= 1e0

    end

    @testset "app algorithm with AC power flow" begin
        data_area = DPM.run_dopf_app(data_WB5, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_14, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_30, "AC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_118, "AC", nlp_solver, 10000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "AC", nlp_solver)
        @test dist_error <= 1e-1

    end

    @testset "app algorithm with SOC relaxation of power flow" begin
        data_area = DPM.run_dopf_app(data_WB5, "SOCP", nlp_solver, 100, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_14, "SOCP", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_30, "SOCP", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e-0

        data_area = DPM.run_dopf_app(data_118, "SOCP", nlp_solver, 100000, tol = 1e-2, verbose = false)
        dist_error = DPM.compare_solution(data_118, data_area, "SOCP", nlp_solver)
        @test dist_error <= 1e0

    end

    @testset "app algorithm with QC relaxation of power flow" begin
        data_area = DPM.run_dopf_app(data_WB5, "QC", nlp_solver, 100, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_WB5, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_14, "QC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_14, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-2

        data_area = DPM.run_dopf_app(data_30, "QC", nlp_solver, 1000, tol = 1e-4, verbose = false)
        dist_error = DPM.compare_solution(data_30, data_area, "QC", nlp_solver)
        @test dist_error <= 1e-0
    end

end
