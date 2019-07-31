struct FeatureMatrix
    matrix::Array{Float64, 2}
    slots::Array{Any, 1}
end

"""
    buildRowArray(matrix, row)::Array{Int64,2}

Builds a `1xsize(matrix)[1]` row matrix with a `1` in position `row` and `0`
elsewhere.

Examples:
```jldoctest
julia> buildRowArray([1 2 3; 4 5 6; 7 8 9], 2)
1×3 Array{Int64,2}:
 0  1  0
```
"""
function buildRowArray(matrix, row)::Array{Int64,2}
    local rowArray = [0 for x in 1:size(matrix)[1]]
    rowArray[row] = 1
    return hcat(rowArray...)
end

"""
    slotFunc(matrix, row, slotFuncName)

Add a function named `slotFuncName` to access a property of a slot from row
`row` in feature matrix `matrix`.
"""
macro slotFunc(matrix, row, slotFuncName)
    quote
        local mat = $(esc(matrix)).matrix
        local slots = $(esc(matrix)).slots
        $(esc(slotFuncName))(i)::GenericAffExpr = (buildRowArray(mat, $(esc(row))) * mat * slots[i])[1]
        $(esc(slotFuncName))()::GenericAffExpr = (buildRowArray(mat, $(esc(row))) * mat * slots[1])[1]
    end
end

function constructMotorFeatureMatrix(motors, gearRatios)
    return hcat([[
    motor.τStall / ratio
    motor.ωFree * ratio
    motor.price
    motor.mass
    ratio
    ] for motor in motors,
        ratio in gearRatios
    ]...)
end

function constructLinkFeatureMatrix(limb::Limb, rangeLength::Int64)
    return hcat([[
    link1
    link2
    link3
    log(link1)
    log(link2)
    log(link3)
    log(link1 + link2 + link3)
    log(link2 + link3)
    ] for link1 in range(limb.minLinks[1].dhParam.r, stop = limb.maxLinks[1].dhParam.r, length = rangeLength),
        link2 in range(limb.minLinks[2].dhParam.r, stop = limb.maxLinks[2].dhParam.r, length = rangeLength),
        link3 in range(limb.minLinks[3].dhParam.r, stop = limb.maxLinks[3].dhParam.r, length = rangeLength)
    ]...)
end

function constructMotorAndLinkFeatureMatrix(motors, gearRatios, limb::Limb, rangeLength::Int64)
    return hcat([[
    motor.τStall / ratio
    motor.ωFree * ratio
    motor.price
    motor.mass
    ratio
    log(motor.ωFree * ratio)
    link1
    link2
    link3
    log(link1)
    log(link2)
    log(link3)
    log(link1 + link2 + link3)
    log(link2 + link3)
    motor.mass * link1
    motor.mass * link2
    getTorque(limb, link1, link2, link3, 2)
    getTorque(limb, link1, link2, link3, 3)
    ] for motor in motors,
        ratio in gearRatios,
        link1 in range(limb.minLinks[1].dhParam.r, stop = limb.maxLinks[1].dhParam.r, length = rangeLength),
        link2 in range(limb.minLinks[2].dhParam.r, stop = limb.maxLinks[2].dhParam.r, length = rangeLength),
        link3 in range(limb.minLinks[3].dhParam.r, stop = limb.maxLinks[3].dhParam.r, length = rangeLength)
    ]...)
end

function getTorque(limb::Limb, link1Length, link2Length, link3Length, i)
    torques = calculateJointTorques(limb, link1Length, link2Length, link3Length)
    if typeof(torques) <: Array
        return torques[i]
    else
        return torques
    end
end

function calculateJointTorques(limb::Limb, link1Length, link2Length, link3Length)
    θ3 = π/2 - atan(limb.targetZ/1000, limb.targetX/1000)
    h = sqrt((limb.targetX/1000)^2 + ((limb.targetZ/1000) - link1Length)^2)

    if h > link2Length + link3Length
        return 1e+8
    end

    θ1asinArg = (link1Length * sin(θ3)) / h
    if abs(θ1asinArg) > 1
        return 1e+8
    end
    θ1 = asin(θ1asinArg)

    βacosArg = (-h^2 + link2Length^2 + link3Length^2) / (2 * link2Length * link3Length)
    if abs(βacosArg) > 1
        return 1e+8
    end
    β = acos(βacosArg)

    θ4asinArg = (link2Length * sin(β)) / h
    if abs(θ4asinArg) > 1
        return 1e+8
    end
    θ4 = asin(θ4asinArg)

    α = (π - β - θ4) + (π - θ1 - θ3)
    J = [0 (link3Length * cos(α + β) - link2Length * sin(α)) (link3Length * cos(α + β));
         (link1Length + link3Length * sin(α + β) + link2Length * cos(α)) 0 0;
         0 (-link3Length * sin(α + β) - link2Length * cos(α)) (-link3Length * sin(α + β))]
    return transpose(J) * [0, 0, limb.tipForceClosePos]
end
