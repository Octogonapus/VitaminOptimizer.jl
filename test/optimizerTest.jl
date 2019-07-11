using Test, VitaminOptimizer

@testset "loadAndOptimize!" begin
    paretoSolutions = loadAndOptimize!(
        makeGLPKModel(),
        "testConstraints1.json",
        "HephaestusArmLimbOne",
        "testMotorOptions.json"
    )

    @test size(paretoSolutions) == (20,3)

    @test length(collect(Set([vec(paretoSolutions[i,:]) for i in 1:size(paretoSolutions)[1]]))) == 20

    @test isequal(
        # Take the first solution. Don't test for gear ratio because it can vary.
        [x[1] for x in paretoSolutions[1,:]],
        [
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12),
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12),
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
        ]
    )

    @testset "loadAndOptimzeAtParetoFrontier!" begin
        fullyOptimalSolution = loadAndOptimzeAtParetoFrontier!(
            makeGLPKModel(),
            "testConstraints1.json",
            "HephaestusArmLimbOne",
            "testMotorOptions.json"
        )

        # Check the motors in the solution on the Pareto frontier
        @test isequal(
            vcat([x[1] for x in fullyOptimalSolution]...),
            [
                VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12),
                VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12),
                VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
            ]
        )

        # Check the sum of the gear ratios in the fully optimal solution is equal to
        # the maximum of the sum of gear ratios in all solutions on the Paret frontier
        @test isequal(
            sum([x[2] for x in fullyOptimalSolution]),
            maximum([sum(x[2] for x in paretoSolutions[i,:]) for i in 1:size(paretoSolutions)[1]])
        )
    end
end
