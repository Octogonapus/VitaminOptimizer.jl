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

const gravity = 9.80665

"""
	findOptimalChoices(featureMatrix, motors, links)

Find the optimal motors, gear ratios, and link lengths.
"""
findOptimalChoices(featureMatrix::FeatureMatrix, motors, links) =
	hcat([(motors[findMotorIndex(col, motors)], col[5], links[1], links[2], links[3])
		for col in optimalColumns(featureMatrix)]...)

"""
	optimizeAtParetoFrontier(model::Model, objectiveFunction,
		featureMatrix::FeatureMatrix, links, motors, filename::String)

Further optimize the `model` at the Pareto frontier defined by the current objective value.
"""
function optimizeAtParetoFrontier(model::Model, objectiveFunction,
	featureMatrix::FeatureMatrix, links, motors, filename::String)
	stayAtParetoFrontier(model, objectiveFunction)

	@objective(model, Max, sum(links))

	optimize!(model)

	if failedToOptimize(model)
		error("Failed to opimize the model at the Pareto frontier. Ending with objective value: " *
			string(objective_value(model)))
	else
		solution = findOptimalChoices(featureMatrix, motors, value.(links))
		printOptimizationResult!(model, solution, value.(links))
		saveOptimizationResult(model, solution, value.(links), filename)
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
	linkRangeLength = 5
	numFlCols = linkRangeLength^length(limb.maxLinks)

	@variable(model, slot1[1:numFmCols], Bin)
	@variable(model, slot2[1:numFmCols], Bin)
	@variable(model, slot3[1:numFmCols], Bin)
	slots = [slot1, slot2, slot3]
    Fm = FeatureMatrix(constructMotorFeatureMatrix(motors, gearRatios), slots)

	@slotFunc(Fm, 1, motorSlotτ)
	@slotFunc(Fm, 2, motorSlotω)
	@slotFunc(Fm, 3, motorSlotPrice)
	@slotFunc(Fm, 4, motorSlotMass)
	@slotFunc(Fm, 5, motorSlotGearRatio)

	@constraint(model, sum(slot1) == 1)
	@constraint(model, sum(slot2) == 1)
	@constraint(model, sum(slot3) == 1)

	@variable(model, limb.minLinks[1].dhParam.r <= link1 <= limb.maxLinks[1].dhParam.r, Int)
	@variable(model, limb.minLinks[2].dhParam.r <= link2 <= limb.maxLinks[2].dhParam.r, Int)
	@variable(model, limb.minLinks[3].dhParam.r <= link3 <= limb.maxLinks[3].dhParam.r, Int)

	@expression(model, τ1Required, limb.tipForce * (link1 / 1000 + link2 / 1000 + link3 / 1000) +
								   gravity * (motorSlotMass(2) * link1 / 1000 +
								   motorSlotMass(3) * (link1 / 1000) + motorSlotMass(3) * (link2 / 1000)))
	@constraint(model, eq3, motorSlotτ(1) .>= τ1Required)

	@expression(model, τ2Required, limb.tipForce * (link2 / 1000 + link3 / 1000) +
								   gravity * motorSlotMass(3) * (link2 / 1000))
	@constraint(model, eq4, motorSlotτ(2) .>= τ2Required)

	@expression(model, τ3Required, limb.tipForce * (link3 / 1000))
	@constraint(model, eq5, motorSlotτ(3) .>= τ3Required)

	@constraint(model, eq6, motorSlotω(1) * (link1 / 1000 + link2 / 1000 + link3 / 1000) .>= limb.tipVelocity)

	@constraint(model, eq7, motorSlotω(2) * (link2 / 1000 + link3 / 1000) .>= limb.tipVelocity)

	@constraint(model, eq8, motorSlotω(3) * (link3 / 1000) .>= limb.tipVelocity)

	@constraint(model, link1 + link2 + link3 == 400)

	objectiveFunction = sum(i -> motorSlotPrice(i), 1:length(slots))
	@objective(model, Min, objectiveFunction)

	# Run the first optimization pass.
	optimize!(model)

	if failedToOptimize(model)
		error("The model was not solved correctly.")
	else
		# Get the first solution and use it to find the other solutions.
		solution = findOptimalChoices(Fm, motors, value.([link1, link2, link3]))

		println("Found solution:")
		printOptimizationResult!(model, solution, value.([link1, link2, link3]))
		saveOptimizationResult(model, solution, value.([link1, link2, link3]), filename)

		return (model, objectiveFunction, solution, Fm, [link1, link2, link3])
	end
end

function printOptimizationResult!(model, optimalMotors, links)
	if termination_status(model) == MOI.TIME_LIMIT
		println("-------------------------------------------------------")
		println("-------------------SUBOPTIMAL RESULT-------------------")
		println("-------------------------------------------------------")
	end

	println("Optimal objective: ", objective_value(model))
	println("Optimal motors:")
	for (mtr, ratio) in optimalMotors
		println("\t", mtr, ", ratio=", ratio, ", link1=", links[1],
			", link2=", links[2], ", link3=", links[3])
	end
end

function saveOptimizationResult(model, solution, links, filename)
	open(filename, "a") do file
		if termination_status(model) == MOI.TIME_LIMIT
			write(file, "-------------------------------------------------------\n")
			write(file, "-------------------SUBOPTIMAL RESULT-------------------\n")
			write(file, "-------------------------------------------------------\n")
		end

		write(file, "Optimal objective: ", string(objective_value(model)), "\n")
		write(file, "Optimal motors:\n")
		for (mtr, ratio) in solution
			write(file, "\t", string(mtr), ", ratio=", string(ratio), ", link1=", string(links[1]),
				", link2=", string(links[2]), ", link3=", string(links[3]), "\n")
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
	limb, motors, gearRatios = loadProblem(constraintsFile, limbName, motorOptionsFile, 1.0)

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
	limb, motors, gearRatios = loadProblem(constraintsFile, limbName, motorOptionsFile, 1.0)

	println("Optimizing initial model.")
	model, objectiveFunction, solution, featureMatrix, links =
		buildAndOptimizeModel!(model, limb, motors, gearRatios, resultsFile)

	println("Optimizing at Pareto frontier.")
	solution = optimizeAtParetoFrontier(model, objectiveFunction, featureMatrix, links, motors, resultsFile)

	return solution
end

solution = loadAndOptimzeAtParetoFrontier!(
    #makeGLPKModel(),
    makeGurobiModel(1),
    "res/constraints2.json",
    "HephaestusArmLimbOne",
    "res/motorOptions.json",
    "optimizationTestResults_miqp_loadAndOptimzeAtParetoFrontier.txt"
)
