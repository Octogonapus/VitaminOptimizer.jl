import JSON

struct Motor
    name::String
    τStall::Float64 # Units Nm
    ωFree::Float64  # Units rad/s
    price::Float64  # Units usd
    mass::Float64   # Units kg
end

function parseMotorOptions!(fileName::String)::Array{Motor, 1}
    return parseMotorOptions(JSON.parsefile(fileName))
end

function parseMotorOptions(json::Dict{String, Any})::Array{Motor, 1}
    motorClasses = json["motors"]
    entries = collect(motorClasses)
    motors = collect(Iterators.flatten(
        [parseMotorClass(json, className, members) for
            (className, members::Array{String, 1})=entries]))
    return motors
end

function parseMotorClass(
    json::Dict{String, Any},
    className::String,
    members::Array{String, 1})
    return [parseMotor(className * "-" * member,
            json["data"][className * "-" * member]) for member=members]
end

function parseMotor(name::String, motorData::Dict{String, Any})
    return Motor(
        name,
        motorData["MaxTorqueNewtonmeters"],
        motorData["MaxFreeSpeedRadPerSec"],
        motorData["price"],
        motorData["massKg"])
end
