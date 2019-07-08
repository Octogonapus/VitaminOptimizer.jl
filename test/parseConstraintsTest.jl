using Test

include("../src/parseConstraints.jl")

@testset "parseConstraints" begin
        @test parseConstraints!("res/testConstraints1.json", ["HephaestusArmLimbOne"]) == Dict(
                "max" => [DhParam(0.135, 0, 0, -90), DhParam(0, 0, 0.175, 0), DhParam(0, 90, 0.16928, 0)],
                "min" => [DhParam(0.135, 0, 0, -90), DhParam(0, 0, 0.175, 0), DhParam(0, 90, 0.16928, 0)])

        @test parseConstraints!("res/testConstraints2.json", ["HephaestusArmLimbOne"]) == Dict(
                "max" => [DhParam(0.135, 0, 0.1, -90), DhParam(0, 0, 0.2, 0), DhParam(0, 90, 0.2, 0)],
                "min" => [DhParam(0.135, 0, 0, -90), DhParam(0, 0, 0.175, 0), DhParam(0, 90, 0.16928, 0)])
end
