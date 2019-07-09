using Test, VitaminOptimizer

@testset "loadAndOptimize!" begin
    @test isequal(
        # Don't test for gear ratio because it can vary right now
        [mtr[1] for mtr in loadAndOptimize!(
            makeGLPKModel(),
            "res/testConstraints1.json",
            "HephaestusArmLimbOne",
            "res/testMotorOptions.json"
        )],
        [
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
            VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12)
        ]
    )
end
