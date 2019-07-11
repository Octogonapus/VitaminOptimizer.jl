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
    println("===========buildRowArray================")
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
        function $(esc(slotFuncName))(i)::GenericAffExpr
            return (buildRowArray(mat, $(esc(row))) * mat * slots[i])[1]
        end
    end
end

function constructMotorFeatureMatrix(motors, gearRatios)
    return hcat([
    [motor.τStall / ratio
    motor.ωFree * ratio
    motor.price
    motor.mass
    (motor.ωFree * ratio) / (motor.τStall / ratio)
    ratio] for motor in motors,
        ratio in gearRatios
   ]...)
end

# function constructMotorFeatureMatrix(motors, gearRatios, limb::Limb, linkRangeLength::Int64)
#     return hcat([
#     [motor.τStall / ratio
#     motor.ωFree * ratio
#     motor.price
#     motor.mass
#     (motor.ωFree * ratio) / (motor.τStall / ratio)
#     ratio
#     link1
#     link2
#     link3
#     motor.mass * link1
#     motor.mass * link2] for motor in motors,
#         ratio in gearRatios,
#         link1 in range(limb.minLinks[1].dhParam.r, stop = limb.maxLinks[1].dhParam.r, length = linkRangeLength),
#         link2 in range(limb.minLinks[2].dhParam.r, stop = limb.maxLinks[2].dhParam.r, length = linkRangeLength),
#         link3 in range(limb.minLinks[3].dhParam.r, stop = limb.maxLinks[3].dhParam.r, length = linkRangeLength)
#    ]...)
# end
