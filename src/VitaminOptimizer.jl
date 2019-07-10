module VitaminOptimizer

import JSON, GLPK
using JuMP, LinearAlgebra

include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

const gravity = 9.80665

"""
	optimalIndices(slots)

Find the indices of the chosen values of `slots` (the indices where `slots[i] == 1`).
"""
optimalIndices(slots) = [findfirst(isequal(1), value.(slot)) for slot in slots]

"""
	optimalColumns(featureMatrix, slots)

Get the columns of the `featureMatrix` corresponding to the chosen values for the `slots`.
"""
optimalColumns(featureMatrix, slots) = [featureMatrix[:, i] for i in optimalIndices(slots)]

"""
	findMotorIndex(featureMatrixColumn, motors)

Find the index of the motor in the `motors` array by searching for a motor with τStall, ωFree, price,
and mass equal to those in the `featureMatrixColumn` (after un-applying the gear ratio).
"""
findMotorIndex(featureMatrixColumn, motors) = findfirst(
	# Approximate equality on τStall and ωFree because we are un-applying the gear ratio.
	x::Motor -> x.τStall ≈ featureMatrixColumn[1] * featureMatrixColumn[6] &&
		x.ωFree ≈ featureMatrixColumn[2] / featureMatrixColumn[6] &&
		x.price == featureMatrixColumn[3] &&
		x.mass == featureMatrixColumn[4],
	motors)

"""
	findOptimalMotors(featureMatrix, allSlots, motors)

Find the optimal motors after the most recent optimization.
"""
findOptimalMotors(featureMatrix, allSlots, motors) = reshape(
	[(motors[findMotorIndex(col, motors)], col[6]) for col in optimalColumns(featureMatrix, allSlots)],
	1, length(allSlots))

"""
	findAllSolutions(model, optimalObjectiveValue, featureMatrix, allSlots, motors)

Iteratively optimize the `model` to find all solutions at a given `optimalObjectiveValue` by adding a
constraint to disallow the most recent combination of slot values. Returns an array of all solutions.
"""
function findAllSolutions(model, optimalObjectiveValue, featureMatrix, allSlots, motors)
	# Disallow the current solution by disallowing the combination of the current slot1, slot2, and slot3
	# values.
	@constraint(model, sum(x -> x[1][x[2]], zip(allSlots, optimalIndices(allSlots))) <= length(allSlots) - 1)

	# Optimize again to find a different solution.
	optimize!(model)

	# If we are no longer at the given optimum or if the model failed to opimize, stop.
	if objective_value(model) != optimalObjectiveValue || failedToOptimize(model)
		return []
	else
		# Record the current solution.
		solution = findOptimalMotors(featureMatrix, allSlots, motors)

		# Keep finding more solutions.
		otherSolutions = findAllSolutions(model, optimalObjectiveValue, featureMatrix, allSlots, motors)

		# Add the current solution to the end of the other solutions.
		if otherSolutions == []
			return solution
		else
			return vcat(solution, otherSolutions)
		end
	end
end

"""
	failedToOptimize(model)

Check if the `model` failed to optimize.
"""
failedToOptimize(model) = !(termination_status(model) == MOI.OPTIMAL ||
	(termination_status(model) == MOI.TIME_LIMIT && has_values(model)))

function buildAndOptimizeModel!(model, limb, limbConfig, motors, gearRatios)
	F_m = constructMotorFeatureMatrix(motors, gearRatios)
	(numRows, numCols) = size(F_m)

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

	# Equation 3
	@expression(model, τ1Required, limb.tipForce * (limbConfig[1].dhParam.r + limbConfig[2].dhParam.r +
	 							   		limbConfig[3].dhParam.r) +
								   gravity * (slotMass(2) * limbConfig[1].dhParam.r +
								   slotMass(3) * (limbConfig[1].dhParam.r + limbConfig[2].dhParam.r)))
	@constraint(model, eq3, slotτ(1) .>= τ1Required)

	# Equation 4
	@expression(model, τ2Required, limb.tipForce * (limbConfig[2].dhParam.r + limbConfig[3].dhParam.r) +
								   slotMass(3) * gravity * limbConfig[2].dhParam.r)
	@constraint(model, eq4, slotτ(2) .>= τ2Required)

	# Equation 5
	@expression(model, τ3Required, limb.tipForce * limbConfig[3].dhParam.r)
	@constraint(model, eq5, slotτ(3) .>= τ3Required)

	# Equation 6
	@expression(model, ω1Required, limb.tipVelocity / (limbConfig[1].dhParam.r + limbConfig[2].dhParam.r +
	 									limbConfig[3].dhParam.r))
	@constraint(model, eq6, slotω(1) .>= ω1Required)

	# Equation 7
	@expression(model, ω2Required, limb.tipVelocity / (limbConfig[2].dhParam.r + limbConfig[3].dhParam.r))
	@constraint(model, eq7, slotω(2) .>= ω2Required)

	# Equation 8
	@expression(model, ω3Required, limb.tipVelocity / limbConfig[3].dhParam.r)
	@constraint(model, eq8, slotω(3) .>= ω3Required)

	@objective(model, Min, sum(x -> slotPrice(x), collect(1:numSlots)))

	# Run the first optimization pass.
	optimize!(model)

	if failedToOptimize(model)
		error("The model was not solved correctly.")
	else
		# Get the first solution and use it to find the other solutions.
		solution = findOptimalMotors(F_m, allSlots, motors)
		otherSolutions = findAllSolutions(model, objective_value(model), F_m, allSlots, motors)
		return (model, vcat(solution, otherSolutions))
	end
end

function makeGLPKModel()::Model
	return Model(with_optimizer(GLPK.Optimizer))
end

function printOptimizationResult!(model, optimalMotors)
	println("Optimal objective: ", objective_value(model))
	println("Optimal motors:")
	for (mtr, ratio) in optimalMotors
		println("\t", mtr, ", ratio=", ratio)
	end
end

"""
	loadAndOptimize!(model, constraints, armName, motorOptions)

Load the constraints from `constraintsFile` and the motor options from
`motorOptionsFile`. Select limb `limbName` from the constraints. Uses the
optimizer in the `model`. Prints and returns the optimal motors.

# Examples
```jldoctest
julia> loadAndOptimize!("res/constraints1.json", "HephaestusArmLimbOne", "res/motorOptions.json")
Optimal objective: 38.849999999999994
Optimal motors:
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.021
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.048
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.077
```
"""
function loadAndOptimize!(model::Model,
		constraintsFile::String,
		limbName::String,
		motorOptionsFile::String)
	limb = parseConstraints!(constraintsFile, [limbName])[1]
	limbConfig = limb.minLinks

	motors = parseMotorOptions!(motorOptionsFile)

	# TODO: Put available gear ratios in the constraints file
	ratios = collect(range(1, step=2, length=30))
	gearRatios = Set(hcat(ratios, 1 ./ ratios))

	(model, allSolutions) = buildAndOptimizeModel!(model, limb, limbConfig, motors, gearRatios)

	if termination_status(model) == MOI.TIME_LIMIT
		println("-------------------------------------------------------")
		println("-------------------SUBOPTIMAL RESULT-------------------")
		println("-------------------------------------------------------")
	end

	return allSolutions
end

export loadAndOptimize!
export makeGLPKModel

end # module VitaminOptimizer
