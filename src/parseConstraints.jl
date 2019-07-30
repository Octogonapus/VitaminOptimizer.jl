import JSON

struct DhParam
    d::Float64 # Units m
    θ::Float64 # Units deg
    r::Float64 # Units m
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
    targetX::Float64     # Units m
    targetZ::Float64     # Units m
end

function parseConstraints!(fileName::String, limbNames::Array{String, 1}, linearConversion::Float64)::Array{Limb, 1}
    return parseConstraints(JSON.parsefile(fileName), limbNames, linearConversion)
end

function parseConstraints(json::Dict{String, Any}, limbNames::Array{String, 1}, linearConversion::Float64)::Array{Limb, 1}
    return [parseLimb(json, name, linearConversion) for name in limbNames]
end

function parseLimb(json::Dict{String, Any}, limbName::String, linearConversion::Float64)
    tipVelocity = json["requiredTipVelocityMeterPerSec"]
    tipForce = json["requiredTipForceNewtons"]
    limb = json[limbName]

    (maxConfig, minConfig) = parseLimbConfig(limb, linearConversion)
    @assert length(maxConfig) == length(minConfig)

    return Limb(
        limbName,
        maxConfig,
        minConfig,
        tipVelocity,
        tipForce,
        0.05,
        0.05
    )
end

function parseLimbConfig(limb::Dict{String, Any}, linearConversion::Float64)
    minConfig = parseLinks(limb["min"], linearConversion)
    maxConfig = parseLinks(limb["max"], linearConversion)
    return (maxConfig, minConfig)
end

function parseLinks(config::Dict{String, Any}, linearConversion::Float64)::Tuple{Vararg{Link}}
    links = sort(collect(values(config)), by=x -> x["index"])
    return tuple(map(x -> Link(parseDhParam(x, linearConversion)), links)...)
end

function parseDhParam(link::Dict{String, Any}, linearConversion::Float64)::DhParam
    return DhParam(link["dh-D"]/linearConversion, link["dh-Theta"], link["dh-A"]/linearConversion, link["dh-Alpha"])
end
