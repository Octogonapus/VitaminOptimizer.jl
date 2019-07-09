import JSON

struct DhParam
    d::Float64     # Units m
    θ::Float64 # Units deg
    r::Float64     # Units m
    α::Float64 # Units deg
end

function Base.show(io::IO, dhParam::DhParam)
    print(
        io,
        "DhParam(d=",
        dhParam.d,
        ", θ=",
        dhParam.θ,
        ", r=",
        dhParam.r,
        ", α=",
        dhParam.α,
        ")"
    )
end

struct Link
    dhParam::DhParam
end

struct Limb
    name::String
    maxLinks::Tuple{Vararg{Link}}
    minLinks::Tuple{Vararg{Link}}
    tipVelocity::Float64 # Units m/s
    tipForce::Float64    # Units N
end

function parseConstraints!(fileName::String, limbNames::Array{String, 1})::Array{Limb, 1}
    return parseConstraints(JSON.parsefile(fileName), limbNames)
end

function parseConstraints(json::Dict{String, Any}, limbNames::Array{String, 1})::Array{Limb, 1}
    return [parseLimb(json, name) for name in limbNames]
end

function parseLimb(json::Dict{String, Any}, limbName::String)
    tipVelocity = json["requiredTipVelocityMeterPerSec"]
    tipForce = json["requiredTipForceNewtons"]
    limb = json[limbName]
    (maxConfig, minConfig) = parseLimbConfig(limb)
    return Limb(
        limbName,
        maxConfig,
        minConfig,
        tipVelocity,
        tipForce
    )
end

function parseLimbConfig(limb::Dict{String, Any})
    minConfig = parseLinks(limb["min"])
    maxConfig = parseLinks(limb["max"])
    return (maxConfig, minConfig)
end

function parseLinks(config::Dict{String, Any})::Tuple{Vararg{Link}}
    links = sort(collect(values(config)), by=x -> x["index"])
    return tuple(map(x -> Link(parseDhParam(x)), links)...)
end

function parseDhParam(link::Dict{String, Any})::DhParam
    return DhParam(link["dh-D"] / 1000, link["dh-Theta"], link["dh-A"] / 1000, link["dh-Alpha"])
end
