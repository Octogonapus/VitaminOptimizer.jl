using Test

include("../src/parseMotorOptions.jl")

@testset "parseMotorOptions" begin
    @test parseMotorOptions!("res/motorOptions.json") == [
        Motor(0.098, 139.626, 0.12, 12.95),
        Motor(0.3108, 87.2665, 0.365, 19.95),
        Motor(1.4123, 20.943, 0.203, 38.95)]
end
