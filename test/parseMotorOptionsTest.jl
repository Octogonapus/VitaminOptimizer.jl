using Test

include("../src/parseMotorOptions.jl")

@testset "parseMotorOptions" begin
    @test parseMotorOptions!("testMotorOptions.json") == [
        Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12),
        Motor("stepperMotor-GenericNEMA17", 0.3108, 87.2665, 19.95, 0.365),
        Motor("roundMotor-WPI-gb37y3530-50en", 1.4123, 20.943, 38.95, 0.203)]
end
