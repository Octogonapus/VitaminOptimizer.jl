import JSON

struct DhParam
    d::Float64     # Units m
    theta::Float64 # Units deg
    r::Float64     # Units m
    alpha::Float64 # Units deg
end

function parseConstraints!(fileName::String, limbNames::Array{String, 1})
    return parseConstraints(JSON.parsefile(fileName), limbNames)
end

function parseConstraints(json::Dict{String, Any}, limbNames::Array{String, 1})
    tipVelocity = json["requiredTipVelocityMeterPerSec"]
    tipForce = json["requiredTipForceNewtons"]
    limb = json["HephaestusArmLimbOne"]
    return parseLimb(limb)
end

function parseLimb(limb::Dict{String, Any})
    minConfig = parseLimbConfig(limb["min"])
    maxConfig = parseLimbConfig(limb["max"])

    return Dict(
        "max" => maxConfig,
        "min" => minConfig
    )
end

function parseLimbConfig(config::Dict{String, Any})
    links = sort(collect(values(config)), by=x -> x["index"])
    return map(parseDhParam, links)
end

function parseDhParam(link::Dict{String, Any})
    return DhParam(link["dh-D"] / 1000, link["dh-Theta"], link["dh-A"] / 1000, link["dh-Alpha"])
end
