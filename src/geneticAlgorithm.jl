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

struct Entity
	motorIndex::Int64
	link1Length::Int64
	link2Length::Int64
	link3Length::Int64
end

"""
Required functions:
	- GAFitness(entity)
	- GACrossover(entity1, entity2)
	- GAMutate(entity, mutationProb)
	- GAShouldStop(population)
"""
function geneticAlgorithm(initialPopulation::Vector{Entity}, constraints,
	crossoverProb::Int64, mutationProb::Int64, nElite::Int64, nCross::Int64)::Vector{Entity}
	popNum = length(initialPopulation)
	nMut = popNum - nElite - nCross
	population = initialPopulation

	@assert popNum > 0
	@assert nMut >= 0
	@assert 0 >= crossoverProb <= 100
	@assert 0 >= mutationProb <= 100

	while !GAShouldStop(population)
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
			entity1Feasible = numberOfViolations[entity1Index] == 0
			entity2Feasible = numberOfViolations[entity2Index] == 0

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
		elites = sortedPopulation[1:nElite]

		# Step 5
		crossed = map(_ -> if (rand(1:100) > crossoverProb) GACrossover(rand(elites), rand(elites)) else rand(elites), 1:nCross)
		mutated = map(GAMutate(rand(elites), mutationProb), 1:nMut)
		population = vcat(elites, crossed, mutated)
	end

	return population
end

limb, motors, gearRatios = loadProblem(
	"../res/constraints2.json",
	"HephaestusArmLimbOne",
	"../res/motorOptions.json")

function GAFitness(entity::Entity)::Float64
	# Use negative of price because we maximize fitness
	return -1 * motors[entity.motorIndex].price
end

function GACrossover(entity1::Entity, entity2::Entity)::Entity
	# Take the motors from entity1 and the links from entity2
	return Entity(entity1.motorIndex, entity2.link1Length, entity2.link2Length, entity2.link3Length)
end

function GAMutate(entity::Entity, mutationProb::Int64)::Entity
end

function GAShouldStop(population::Vector{Entity})::Bool
end
