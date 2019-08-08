import JSON

include("motor.jl")

function parseMotorOptions!(fileName::String)::Array{Motor, 1}
    return parseMotorOptions(JSON.parsefile(fileName))
end

function parseMotorOptions(json::Array{Any, 1})::Array{Motor, 1}
    return [parseMotor(it) for it=json]
end

function parseMotor(data::Dict{String, Any})
    return Motor(
        data["name"],
        data["τStall"],
        data["ωFree"],
        data["price"],
        data["mass"])
end
