module VitaminOptimizer

import JSON, GLPK
using JuMP, Gurobi, LinearAlgebra

include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

const gravity = 9.80665

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

	optimize!(model)

	function findMotorIndex(featureMatrixColumn)
		return findfirst(
			x::Motor -> x.τStall ≈ featureMatrixColumn[1] * featureMatrixColumn[6] &&
				x.ωFree ≈ featureMatrixColumn[2] / featureMatrixColumn[6] &&
				x.price == featureMatrixColumn[3] &&
				x.mass == featureMatrixColumn[4],
			motors)
	end

	function findOptimalMotors()
		optimalMotorIndices = [findfirst(isequal(1), value.(slot)) for slot in allSlots]
		optimalMotorColumns = [F_m[:, i] for i in optimalMotorIndices]
		optimalMotors = [(motors[findMotorIndex(col)], col[6]) for col in optimalMotorColumns]
		return optimalMotors
	end


	function printOptimizationResult!(optimalMotors)
		println("Optimal objective: ", objective_value(model))
		println("Optimal motors:")
		for (mtr, ratio) in optimalMotors
			println("\t", mtr, ", ratio=", ratio)
		end
	end

	if !(termination_status(model) == MOI.OPTIMAL || (termination_status(model) == MOI.TIME_LIMIT && has_values(model)))
		error("The model was not solved correctly.")
	else
		optimalMotors = findOptimalMotors()

		if termination_status(model) == MOI.TIME_LIMIT
			println("-------------------------------------------------------")
			println("-------------------SUBOPTIMAL RESULT-------------------")
			println("-------------------------------------------------------")
		end

		printOptimizationResult!(optimalMotors)
		return (model, optimalMotors)
	end
end

function makeGurobiModel!(maxNumSolutions::Int64)::Model
	env = Gurobi.Env()
	setparam!(env, "PoolSearchMode", 2)
	setparam!(env, "PoolSolutions", 10)
	return Model(with_optimizer(Gurobi.Optimizer, env, Presolve=1))
end

function makeGLPKModel()::Model
	return Model(with_optimizer(GLPK.Optimizer))
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
		motorOptionsFile::String
	)
	limb = parseConstraints!(constraintsFile, [limbName])[1]
	limbConfig = limb.minLinks

	motors = parseMotorOptions!(motorOptionsFile)
	ratios = collect(range(1, step=2, length=30))
	gearRatios = Set(hcat(ratios, 1 ./ ratios))

	(model, optimalMotors) = buildAndOptimizeModel!(model, limb, limbConfig, motors, gearRatios)
	return optimalMotors
end

export loadAndOptimize!
export makeGurobiModel!
export makeGLPKModel

end # module VitaminOptimizer
