using Test

include("../src/parseMotorOptions.jl")
include("../src/featureMatrix.jl")

function featureMatrixToSet(fm)
    return Set([Set(fm[:,n]) for n=1:size(fm)[2]])
end

@testset "constructMotorFeatureMatrix" begin
    @test featureMatrixToSet(
        constructMotorFeatureMatrix(
            parseMotorOptions!("res/testMotorOptions.json"),
            [2 1 1/3])
        ) == featureMatrixToSet(hcat(
            [0.098 / 2
             139.626 * 2
             12.95
             0.12
             (139.626 * 2) / (0.098 / 2)],

            [0.098
             139.626
             12.95
             0.12
             139.626 / 0.098],

            [0.098 / (1/3)
             139.626 * (1/3)
             12.95
             0.12
             (139.626 * (1/3)) / (0.098 / (1/3))],

            [0.3108 / 2
             87.2665 * 2
             19.95
             0.365
             (87.2665 * 2) / (0.3108 / 2)],

            [0.3108
             87.2665
             19.95
             0.365
             87.2665 / 0.3108],

            [0.3108 / (1/3)
             87.2665 * (1/3)
             19.95
             0.365
             (87.2665 * (1/3)) / (0.3108 / (1/3))],

            [1.4123 / 2
             20.943 * 2
             38.95
             0.203
             (20.943 * 2) / (1.4123 / 2)],

            [1.4123
             20.943
             38.95
             0.203
             20.943 / 1.4123],

            [1.4123 / (1/3)
             20.943 * (1/3)
             38.95
             0.203
             (20.943 * (1/3)) / (1.4123 / (1/3))]))
end
