import JSON
using Distributions

include("motor.jl")

function writeMotorsToFile!(fileName::String, motors::Array{Motor, 1})
    open(fileName, "w") do outFile
        JSON.print(outFile, motors)
    end
end

function generateRandomMotor()::Motor
    return Motor(
        string(rand(Uniform(1, 1e+5))),
        rand(Uniform(0.01, 40)),
        rand(Uniform(5, 35)),
        rand(Uniform(0.01, 40)),
        rand(Uniform(0.01, 0.2))
    )
end

writeMotorsToFile!(
    "res/random_motor_options_50.json",
    [generateRandomMotor() for i=1:50]
)
