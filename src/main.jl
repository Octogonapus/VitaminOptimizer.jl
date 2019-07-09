import JSON, GLPK
using JuMP
# using Gurobi

include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

limbConstraints = parseConstraints!("res/constraints1.json", ["HephaestusArmLimbOne"])
limbConfig = limbConstraints["min"]

motors = parseMotorOptions!("res/motorOptions.json")
gearRatios = [7, 5, 3, 1, 1/3, 1/5, 1/7]
F_m = constructMotorFeatureMatrix(motors, gearRatios)

const gravity = 9.80665

# env = Gurobi.Env()
# setparam!(env, "LogFile",
# 		  "/home/salmon/Documents/auto-configured-vitamins-optimizer/main.log")
#
# model = Model(with_optimizer(Gurobi.Optimizer, Presolve=1))
model = Model(with_optimizer(GLPK.Optimizer))

(numRows, numCols) = size(F_m)

τRow = [1 0 0 0 0]
ωRow = [0 1 0 0 0]
priceRow = [0 0 1 0 0]
massRow = [0 0 0 1 0]
omegaFuncRow = [0 0 0 0 1]

# Each slot is a binary vector with a 1 that picks which motor to use.
@variable(model, slot1[1:numCols], Bin)
@constraint(model, slot1Unique, sum(slot1) == 1)

@variable(model, slot2[1:numCols], Bin)
@constraint(model, slot2Unique, sum(slot2) == 1)

@variable(model, slot3[1:numCols], Bin)
@constraint(model, slot3Unique, sum(slot3) == 1)

allSlots = [slot1, slot2, slot3]

"""
	@ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i::Int64)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i::Int64) =
	((τRow * F_m * allSlots[i]) - τ) * (omegaFuncRow * F_m * allSlots[i])

"""
	@ω(τ::Float64, i::Int64)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Float64, i::Int64) =
	((τRow * F_m * allSlots[i]) .- τ) * (omegaFuncRow * F_m * allSlots[i])

# Equation 3
# limbConfig[1].alpha is 0 which makes this boring
@expression(model, τ1Required, tipForce * (limbConfig[1].alpha + limbConfig[2].alpha + limbConfig[3].alpha) +
							   gravity * (massRow * F_m * slot2 * limbConfig[1].alpha +
							   massRow * F_m * slot3 * (limbConfig[1].alpha + limbConfig[2].alpha)))
@constraint(model, eq3, τRow * F_m * slot1 .>= τ1Required)

# Equation 4
@expression(model, τ2Required, tipForce * (limbConfig[2].alpha + limbConfig[3].alpha) +
							   massRow * F_m * slot3 * gravity * limbConfig[2].alpha)
@constraint(model, eq4, τRow * F_m * slot2 .>= τ2Required)

# Equation 5
@expression(model, τ3Required, tipForce * limbConfig[3].alpha)
@constraint(model, eq5, τRow * F_m * slot3 .>= τ3Required)

# Equation 6
@expression(model, ω1Required, tipVelocity / (limbConfig[1].alpha + limbConfig[2].alpha + limbConfig[3].alpha))
@constraint(model, eq6, ωRow * F_m * slot1 .>= ω1Required)

# Equation 7
@expression(model, ω2Required, tipVelocity / (limbConfig[2].alpha + limbConfig[3].alpha))
@constraint(model, eq7, ωRow * F_m * slot2 .>= ω2Required)

# Equation 8
@expression(model, ω3Required, tipVelocity / limbConfig[3].alpha)
@constraint(model, eq8, ωRow * F_m * slot3 .>= ω3Required)

@objective(model, Min, sum(x -> priceRow * F_m * x, allSlots)[1])

optimize!(model)

println("Optimal objective: ", objective_value(model),
	". slot1 = ", value.(slot1), ", slot2 = ", value.(slot2),
	", slot3 = ", value.(slot3))

optimalMotorIndices = [findfirst(isequal(1), value.(slot)) for slot in allSlots]
optimalMotorNames = [motorNamesWithRatios[i] for i in optimalMotorIndices]
println("Optimal motors: ", optimalMotorNames)
