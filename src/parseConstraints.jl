import JSON

struct DhParam
    d::Float64
    theta::Float64
    r::Float64
    alpha::Float64
end

function parseConstraints!(fileName::String, limbNames::Array{String, 1})
    return parseConstraints(JSON.parsefile(fileName), limbNames)
end

function parseConstraints(json, limbNames::Array{String, 1})
    tipVelocity = json["requiredTipVelocityMeterPerSec"]
    tipForce = json["requiredTipForceNewtons"]
    limb = json["HephaestusArmLimbOne"]
    return parseLimb(limb)
end

function parseLimb(limb)
    minConfig = parseLimbConfig(limb["min"])
    maxConfig = parseLimbConfig(limb["max"])

    return Dict(
        "max" => maxConfig,
        "min" => minConfig
    )
end

function parseLimbConfig(config)
    links = sort(collect(values(config)), by=x -> x["index"])
    return map(parseDhParam, links)
end

function parseDhParam(link)::DhParam
    return DhParam(link["dh-D"] / 1000, link["dh-Theta"], link["dh-A"] / 1000, link["dh-Alpha"])
end
