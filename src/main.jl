import JSON, GLPK
using JuMP, LinearAlgebra
using Gurobi

include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

limb = parseConstraints!("res/constraints1.json", ["HephaestusArmLimbOne"])
limbConfig = limb.minLinks

motors = parseMotorOptions!("res/motorOptions.json")
ratios = collect(range(1, step=2, length=30))
gearRatios = Set(hcat(ratios, 1 ./ ratios))
F_m = constructMotorFeatureMatrix(motors, gearRatios)
(numRows, numCols) = size(F_m)

# model = Model(with_optimizer(Gurobi.Optimizer, Presolve=1))
model = Model(with_optimizer(GLPK.Optimizer))

# Each slot is a binary vector with a 1 that picks which motor to use.
@variable(model, slot1[1:numCols], Bin)
@constraint(model, slot1Unique, sum(slot1) == 1)

@variable(model, slot2[1:numCols], Bin)
@constraint(model, slot2Unique, sum(slot2) == 1)

@variable(model, slot3[1:numCols], Bin)
@constraint(model, slot3Unique, sum(slot3) == 1)

allSlots = [slot1, slot2, slot3]
numSlots = length(allSlots)

featureIdentity = Matrix{Float64}(I, numRows, numRows)
τRow = transpose(featureIdentity[1,:])
slotτ(i) = τRow * F_m * allSlots[i]

ωRow = transpose(featureIdentity[2,:])
slotω(i) = ωRow * F_m * allSlots[i]

priceRow = transpose(featureIdentity[3,:])
slotPrice(i) = priceRow * F_m * allSlots[i]

massRow = transpose(featureIdentity[4,:])
slotMass(i) = massRow * F_m * allSlots[i]

omegaFuncRow = transpose(featureIdentity[5,:])
slotOmegaFunc(i) = omegaFuncRow * F_m * allSlots[i]

gearRatioRow = transpose(featureIdentity[6,:])
slotGearRatioFunc(i) = gearRatioRow * F_m * allSlots[i]

"""
	@ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i) = (slotτ(i) - τ) * slotOmegaFunc(i)

"""
	@ω(τ::Float64, i)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Float64, i) = (slotτ(i) .- τ) * slotOmegaFunc(i)

const gravity = 9.80665

# Equation 3
# limbConfig[1].dhParam.alpha is 0 which makes this boring
@expression(model, τ1Required, limb.tipForce * (limbConfig[1].dhParam.alpha + limbConfig[2].dhParam.alpha +
 							   		limbConfig[3].dhParam.alpha) +
							   gravity * (slotMass(2) * limbConfig[1].dhParam.alpha +
							   slotMass(3) * (limbConfig[1].dhParam.alpha + limbConfig[2].dhParam.alpha)))
@constraint(model, eq3, slotτ(1) .>= τ1Required)

# Equation 4
@expression(model, τ2Required, limb.tipForce * (limbConfig[2].dhParam.alpha + limbConfig[3].dhParam.alpha) +
							   slotMass(3) * gravity * limbConfig[2].dhParam.alpha)
@constraint(model, eq4, slotτ(2) .>= τ2Required)

# Equation 5
@expression(model, τ3Required, limb.tipForce * limbConfig[3].dhParam.alpha)
@constraint(model, eq5, slotτ(3) .>= τ3Required)

# # Equation 6
# @expression(model, ω1Required, limb.tipVelocity / (limbConfig[1].dhParam.alpha + limbConfig[2].dhParam.alpha +
#  									limbConfig[3].dhParam.alpha))
# @constraint(model, eq6, slotω(1) .>= ω1Required)
#
# # Equation 7
# @expression(model, ω2Required, limb.tipVelocity / (limbConfig[2].dhParam.alpha + limbConfig[3].dhParam.alpha))
# @constraint(model, eq7, slotω(2) .>= ω2Required)
#
# # Equation 8
# @expression(model, ω3Required, limb.tipVelocity / limbConfig[3].dhParam.alpha)
# @constraint(model, eq8, slotω(3) .>= ω3Required)

@objective(model, Min, sum(x -> slotPrice(x), collect(1:numSlots)))

optimize!(model)

# println("Optimal objective: ", objective_value(model))
# println("slot1 = ", value.(slot1))
# println("slot2 = ", value.(slot2))
# println("slot3 = ", value.(slot3))
#
# optimalMotorIndices = [findfirst(isequal(1), value.(slot)) for slot in allSlots]
# optimalMotors = [motors[i] for i in optimalMotorIndices]
# println("Optimal motors: ", optimalMotors)

function printOptimizationResult()
	# [println("slot", i, " = ", value.(allSlots[i])) for i in 1:numSlots]
	optimalMotorIndices = [findfirst(isequal(1), value.(slot)) for slot in allSlots]
	optimalMotors = [F_m[:, i] for i in optimalMotorIndices]
	println("Optimal objective: ", objective_value(model))
	println("Optimal motors: ", optimalMotors)
end

if termination_status(model) == MOI.OPTIMAL
	printOptimizationResult()
elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
	println("-------------------------------------------------------")
	println("-------------------SUBOPTIMAL RESULT-------------------")
	println("-------------------------------------------------------")
	printOptimizationResult()
else
    error("The model was not solved correctly.")
end
