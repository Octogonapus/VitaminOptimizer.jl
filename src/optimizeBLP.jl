import JSON, GLPK
using JuMP, LinearAlgebra

include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")
include("featureMatrixUtil.jl")
include("problemUtil.jl")
include("jumpUtil.jl")
include("paretoUtil.jl")
include("gurobiModel.jl")

"""
	findOptimalChoices(featureMatrix, motors)

Find the optimal motors, gear ratios, and link lengths.
"""
findOptimalChoices(featureMatrix::FeatureMatrix, motors) =
	hcat([(motors[findMotorIndex(col, motors)], col[5], col[7], col[8], col[9])
		for col in optimalColumns(featureMatrix)]...)

"""
	optimizeAtParetoFrontier(model::Model, objectiveFunction,
		featureMatrix::FeatureMatrix, motors, filename::String)

Further optimize the `model` at the Pareto frontier defined by the current objective value.
"""
function optimizeAtParetoFrontier(model::Model, objectiveFunction,
	featureMatrix::FeatureMatrix, motors, filename::String)
	stayAtParetoFrontier(model, objectiveFunction)

	@slotFunc(featureMatrix, 7, limbSlotLink1)
	@slotFunc(featureMatrix, 8, limbSlotLink2)
	@slotFunc(featureMatrix, 9, limbSlotLink3)

	@objective(model, Max, limbSlotLink1()*1000 + limbSlotLink2()*1000 + limbSlotLink3()*1000)

	optimize!(model)

	if failedToOptimize(model)
		error("Failed to opimize the model at the Pareto frontier. Ending with objective value: " *
			string(objective_value(model)))
	else
		solution = findOptimalChoices(featureMatrix, motors)
		printOptimizationResult!(model, solution)
		saveOptimizationResult(model, solution, filename)
		return solution
	end
end

"""
	buildAndOptimizeModel!(model, limb, motors, gearRatios)

Add the initial variables and constraints to the `model` using a feature matrix
built from `limb`, and the coproduct of `motors` and `gearRatios`. Optimize the
model to minimize price using the optimizer in the `model`.
"""
function buildAndOptimizeModel!(model::Model, limb::Limb, motors, gearRatios, filename::String)
	limbConfig = limb.minLinks
	numFmCols = length(motors) * length(gearRatios)
	linkRangeLength = 10
	numFlCols = linkRangeLength^length(limb.maxLinks)
	gravity = 9.80665

	@variable(model, slot1[1:numFmCols * numFlCols], Bin)
	@variable(model, slot2[1:numFmCols * numFlCols], Bin)
	@variable(model, slot3[1:numFmCols * numFlCols], Bin)
	slots = [slot1, slot2, slot3]
	Fml = FeatureMatrix(constructMotorAndLinkFeatureMatrix(motors, gearRatios, limb, linkRangeLength), slots)

	@slotFunc(Fml, 1, motorSlotτ)
	@slotFunc(Fml, 2, motorSlotω)
	@slotFunc(Fml, 3, motorSlotPrice)
	@slotFunc(Fml, 4, motorSlotMass)
	@slotFunc(Fml, 5, motorSlotGearRatio)
	@slotFunc(Fml, 6, motorSlotLnω)
	@slotFunc(Fml, 7, limbSlotLink1)
	@slotFunc(Fml, 8, limbSlotLink2)
	@slotFunc(Fml, 9, limbSlotLink3)
	@slotFunc(Fml, 10, limbSlotLnLink1)
	@slotFunc(Fml, 11, limbSlotLnLink2)
	@slotFunc(Fml, 12, limbSlotLnLink3)
	@slotFunc(Fml, 13, limbSlotLnLink123)
	@slotFunc(Fml, 14, limbSlotLnLink23)
	@slotFunc(Fml, 15, massTimesLink1)
	@slotFunc(Fml, 16, massTimesLink2)
    @slotFunc(Fml, 17, closeTorque2)
	@slotFunc(Fml, 17, closeTorque3)

	@constraint(model, sum(slot1) == 1)
	@constraint(model, sum(slot2) == 1)
	@constraint(model, sum(slot3) == 1)
	@constraint(model, limbSlotLink1(1) == limbSlotLink1(2))
	@constraint(model, limbSlotLink1(1) == limbSlotLink1(3))
	@constraint(model, limbSlotLink2(1) == limbSlotLink2(2))
	@constraint(model, limbSlotLink2(1) == limbSlotLink2(3))
	@constraint(model, limbSlotLink3(1) == limbSlotLink3(2))
	@constraint(model, limbSlotLink3(1) == limbSlotLink3(3))

	@expression(model, τ1Required, limb.tipForce * (limbSlotLink1() + limbSlotLink2() + limbSlotLink3()))
	@constraint(model, eq3, motorSlotτ(1) .>= τ1Required)

	@expression(model, τ2Required, limb.tipForce * (limbSlotLink2() + limbSlotLink3()) +
								   gravity * massTimesLink2(3))
	@constraint(model, eq4, motorSlotτ(2) .>= τ2Required)

	@expression(model, τ3Required, limb.tipForce * limbSlotLink3())
	@constraint(model, eq5, motorSlotτ(3) .>= τ3Required)

	@expression(model, ω1Required, log(limb.tipVelocity) - limbSlotLnLink123())
	@constraint(model, eq6, motorSlotLnω(1) .>= ω1Required)

	@expression(model, ω2Required, log(limb.tipVelocity) - limbSlotLnLink23())
	@constraint(model, eq7, motorSlotLnω(2) .>= ω2Required)

	@expression(model, ω3Required, log(limb.tipVelocity) - limbSlotLnLink3())
	@constraint(model, eq8, motorSlotLnω(3) .>= ω3Required)

	@constraint(model, limbSlotLink1() + limbSlotLink2() + limbSlotLink3() == 400 / 1000)

	@constraint(model, motorSlotτ(2) .>= closeTorque2())
	@constraint(model, motorSlotτ(3) .>= closeTorque3())

	objectiveFunction = sum(i -> motorSlotPrice(i), 1:length(slots))
	@objective(model, Min, objectiveFunction)

	# Run the first optimization pass.
	optimize!(model)

	if failedToOptimize(model)
		error("The model was not solved correctly.")
	else
		# Get the first solution and use it to find the other solutions.
		solution = findOptimalChoices(Fml, motors)

		println("Found solution:")
		printOptimizationResult!(model, solution)
		saveOptimizationResult(model, solution, filename)

		return (model, objectiveFunction, solution, Fml)
	end
end

function printOptimizationResult!(model, optimalMotors)
	if termination_status(model) == MOI.TIME_LIMIT
		println("-------------------------------------------------------")
		println("-------------------SUBOPTIMAL RESULT-------------------")
		println("-------------------------------------------------------")
	end

	println("Optimal objective: ", objective_value(model))
	println("Optimal motors:")
	for (mtr, ratio, link1, link2, link3) in optimalMotors
		println("\t", mtr, ", ratio=", ratio, ", link1=", link1*1000,
			", link2=", link2*1000, ", link3=", link3*1000)
	end
end

function saveOptimizationResult(model, solution, filename)
	open(filename, "a") do file
		if termination_status(model) == MOI.TIME_LIMIT
			write(file, "-------------------------------------------------------\n")
			write(file, "-------------------SUBOPTIMAL RESULT-------------------\n")
			write(file, "-------------------------------------------------------\n")
		end

		write(file, "Optimal objective: ", string(objective_value(model)), "\n")
		write(file, "Optimal motors:\n")
		for (mtr, ratio, link1, link2, link3) in solution
			write(file, "\t", string(mtr), ", ratio=", string(ratio), ", link1=", string(link1*1000),
				", link2=", string(link2*1000), ", link3=", string(link3*1000), "\n")
		end
	end
end

"""
	loadAndOptimize!(model::Model, constraintsFile::String, limbName::String,
		motorOptionsFile::String, resultsFile::String)

Load the constraints from `constraintsFile` and the motor options from
`motorOptionsFile`. Select limb `limbName` from the constraints. Uses the
optimizer in the `model`. Returns the Pareto set.

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
function loadAndOptimize!(model::Model, constraintsFile::String, limbName::String,
		motorOptionsFile::String, resultsFile::String)
	limb, motors, gearRatios = loadProblem(constraintsFile, limbName, motorOptionsFile, 1000.0)

	println("Optimizing initial model.")
	model, objectiveFunction, solution, featureMatrix =
		buildAndOptimizeModel!(model, limb, motors, gearRatios, resultsFile)

	println("Exploring Pareto frontier.")
	stayAtParetoFrontier(model, objectiveFunction)
	otherSolutions = exploreParetoFrontier(model, featureMatrix, motors, resultsFile)

	if (isempty(otherSolutions))
		return vcat(solution)
	else
		return vcat(solution, otherSolutions)
	end
end

"""
	loadAndOptimzeAtParetoFrontier!(model::Model, constraintsFile::String,
		limbName::String, motorOptionsFile::String, resultsFile::String)

Load the constraints from `constraintsFile` and the motor options from
`motorOptionsFile`. Select limb `limbName` from the constraints. Uses the
optimizer in the `model`. Optimizes the model once, assumes the objective value from that optimization
is in the Pareto set, and then optimizes again inside the Pareto set. Returns the optimal value from the
second round of optimization.

# Examples
```jldoctest
julia> loadAndOptimize!("res/constraints1.json", "HephaestusArmLimbOne", "res/motorOptions.json")
Optimal objective: 38.849999999999994
Optimal motors:
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.047619
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.047619
    VitaminOptimizer.Motor("stepperMotor-GenericNEMA14", 0.098, 139.626, 12.95, 0.12), ratio=0.111111
```
"""
function loadAndOptimzeAtParetoFrontier!(model::Model, constraintsFile::String,
	limbName::String, motorOptionsFile::String, resultsFile::String)
	limb, motors, gearRatios = loadProblem(
		constraintsFile,
		limbName,
		motorOptionsFile,
		1000.0)

	println("Optimizing initial model.")
	model, objectiveFunction, solution, featureMatrix = buildAndOptimizeModel!(
		model,
		limb,
		motors,
		gearRatios,
		resultsFile)

	println("Optimizing at Pareto frontier.")
	solution = optimizeAtParetoFrontier(
		model,
		objectiveFunction,
		featureMatrix,
		motors,
		resultsFile)

	return solution
end

function timedOptimize(motorOptionsFile::String)
	println("Motor options: ", motorOptionsFile)

	limb, motors, gearRatios = loadProblem(
		"res/constraints2.json",
		"HephaestusArmLimbOne",
		motorOptionsFile,
		1000.0)

	model, objectiveFunction, solution, featureMatrix = @time buildAndOptimizeModel!(
		makeGurobiModel(1),
		limb,
		motors,
		gearRatios,
		"optimizationTestResults_blp_loadAndOptimzeAtParetoFrontier.txt")
end

for it=["res/random_motor_options_10.json",
		"res/random_motor_options_50.json",
		"res/random_motor_options_100.json",
		"res/random_motor_options_500.json",
		"res/random_motor_options_1000.json"]
	timedOptimize(it)
end
