using Test

include("../src/parseConstraints.jl")

@testset "parseConstraints" begin
        @test isequal(
                parseConstraints!("testConstraints1.json", ["HephaestusArmLimbOne"]),
                [Limb(
                        "HephaestusArmLimbOne",
                        tuple([Link(DhParam(0.135, 0, 0, -90)),
                                Link(DhParam(0, 0, 0.175, 0)),
                                Link(DhParam(0, 90, 0.16928, 0))]...),
                        tuple([Link(DhParam(0.135, 0, 0, -90)),
                                Link(DhParam(0, 0, 0.175, 0)),
                                Link(DhParam(0, 90, 0.16928, 0))]...),
                        1,
                        4.90332
                )]
        )

        @test isequal(
                parseConstraints!("testConstraints2.json", ["HephaestusArmLimbOne"]),
                [Limb(
                        "HephaestusArmLimbOne",
                        tuple([Link(DhParam(0.135, 0, 0.1, -90)),
                                Link(DhParam(0, 0, 0.2, 0)),
                                Link(DhParam(0, 90, 0.2, 0))]...),
                        tuple([Link(DhParam(0.135, 0, 0, -90)),
                                Link(DhParam(0, 0, 0.150, 0)),
                                Link(DhParam(0, 90, 0.1, 0))]...),
                        1,
                        4.90332
                )]
        )
end
