using Test, VitaminOptimizer

@testset "loadAndOptimize!" begin
    solutions = loadAndOptimize!(
        makeGLPKModel(),
        "res/testConstraints1.json",
        "HephaestusArmLimbOne",
        "res/testMotorOptions.json"
    )

    for i in 1:size(solutions)[1]
        solution = solutions[i,:]
        println("Solution:")
        for mtr in solution
            println("\t", mtr)
        end
    end

    @test size(solutions) == (20,3)

    @test length(collect(Set([vec(solutions[i,:]) for i in 1:size(solutions)[1]]))) == 20

    @test isequal(
        # Take the first solution. Don't test for gear ratio because it can vary.
        [mtr[1] for mtr in solutions[1,:]],
        [
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
        ]
    )
end
