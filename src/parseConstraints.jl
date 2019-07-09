import JSON

struct DhParam
    d::Float64     # Units m
    theta::Float64 # Units deg
    r::Float64     # Units m
    alpha::Float64 # Units deg
end

struct Link
    dhParam::DhParam
end

struct Limb
    maxLinks::Tuple{Vararg{Link}}
    minLinks::Tuple{Vararg{Link}}
    tipVelocity::Float64 # Units m/s
    tipForce::Float64    # Units N
end

function parseConstraints!(fileName::String, limbNames::Array{String, 1})::Limb
    return parseConstraints(JSON.parsefile(fileName), limbNames)
end

function parseConstraints(json::Dict{String, Any}, limbNames::Array{String, 1})::Limb
    tipVelocity = json["requiredTipVelocityMeterPerSec"]
    tipForce = json["requiredTipForceNewtons"]
    limb = json["HephaestusArmLimbOne"]
    (maxConfig, minConfig) = parseLimb(limb)
    return Limb(
        maxConfig,
        minConfig,
        tipVelocity,
        tipForce
    )
end

function parseLimb(limb::Dict{String, Any})
    minConfig = parseLimbConfig(limb["min"])
    maxConfig = parseLimbConfig(limb["max"])
    return (maxConfig, minConfig)
end

function parseLimbConfig(config::Dict{String, Any})::Tuple{Vararg{Link}}
    links = sort(collect(values(config)), by=x -> x["index"])
    return tuple(map(x -> Link(parseDhParam(x)), links)...)
end

function parseDhParam(link::Dict{String, Any})::DhParam
    return DhParam(link["dh-D"] / 1000, link["dh-Theta"], link["dh-A"] / 1000, link["dh-Alpha"])
end
