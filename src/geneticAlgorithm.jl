include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

const gravity = 9.80665

"""
	loadProblem(constraintsFile::String, limbName::String, motorOptionsFile::String)

Load the constraints from `constraintsFile` and the motor options from
`motorOptionsFile`. Select limb `limbName` from the constraints.
"""
function loadProblem(constraintsFile::String, limbName::String, motorOptionsFile::String)
	limb = parseConstraints!(constraintsFile, [limbName])[1]
	motors = parseMotorOptions!(motorOptionsFile)

	# TODO: Put available gear ratios in the constraints file
	ratios = collect(range(1, step=2, length=10))
	gearRatios = Set(hcat(ratios, 1 ./ ratios))

	return (limb, motors, gearRatios)
end

"""
Required functions:
	- GAFitness(entity)
	- GAIsFeasible(entity)
	- GACrossover(entity1, entity2)
	- GAMutate(entity)
"""
function geneticAlgorithm(initialPopulation, constraints, crossoverProb::Float64, mutationProb::Float64, Nelite::Int64, Ncross::Int64)
	PopNum = length(initialPopulation)
	Nmut = PopNum - Nelite - Ncross
	population = initialPopulation

	# Step 2
	fitness = map(GAFitness, population)
	constraintValues = map(x -> x[0](x[1]), zip(constraints, population))

	# Equation 13
	constraintViolationValues = sum(x -> max(0, x), constraintValues)

	# Equation 14
	numberOfViolations = sum(x -> x > 0, constraintValues) / length(constraints)

	function customLessThan(entity1, entity2)
		entity1Index = findfirst(isequal(entity1), population)
		entity2Index = findfirst(isequal(entity2), population)
		entity1Feasible = GAIsFeasible(entity1)
		entity2Feasible = GAIsFeasible(entity2)

		if entity1Feasible && entity2Feasible
			# Winner is the one with the highest fitness value
			if fitness[entity1Index] > fitness[entity2Index]
				return entity1
			else
				return entity2
			end
		end

		if xor(entity1Feasible, entity2Feasible)
			# Winner is the feasible one
			if entity1Feasible
				return entity1
			else
				return entity2
			end
		end

		# Both are infeasible
		if numberOfViolations[entity1Index] == numberOfViolations[entity2Index]
			# Winner is the one with the lowest CV
			if constraintViolationValues[entity1Index] < constraintViolationValues[entity2Index]
				return entity1
			else
				return entity2
			end
		else
			# Winner is the one with the lowest NV
			if numberOfViolations[entity1Index] < numberOfViolations[entity2Index]
				return entity1
			else
				return entity2
			end
		end
	end

	# Step 3
	sortedPopulation = sort(population, lt=customLessThan)

	# Step 4
	elites = sortedPopulation[1:Nelite]

	# Step 5
	crossed = map(GACrossover(rand(elites), rand(elites)), 1:Ncross)
	mutated = map(GAMutate(rand(elites)), 1:Nmut)
	population = vcat(elites, crossed, mutated)
end

function GAFitness(entity)
end

function GAIsFeasible(entity)
end

function GACrossover(entity1, entity2)
end

function GAMutate(entity)
end
