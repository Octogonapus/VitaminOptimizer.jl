include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("featureMatrix.jl")

import LinearAlgebra, Distributed

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
	motor1Index::Int64
	motor2Index::Int64
	motor3Index::Int64
	link1Length::Float64
	link2Length::Float64
	link3Length::Float64
	gearRatio1::Float64
	gearRatio2::Float64
	gearRatio3::Float64
end

"""
This implements the algorithm described in:
	Chehouri, A., Younes, R., Perron, J., & Ilinca, A. (2016). A constraint-handling technique
	for genetic algorithms using a violation factor. arXiv preprint arXiv:1610.00976.

Required functions:
	- GAFitness(entity)
	- GACrossover(entity1, entity2)
	- GAMutate(entity, mutationProb)
	- GAShouldStop(population, generationNumber)
"""
function geneticAlgorithm(initialPopulation::Vector{Entity}, constraints,
	mutationProb::Float64, nElite::Int64, nCross::Int64)::Vector{Entity}
	popNum = length(initialPopulation)
	nMut = popNum - nElite - nCross
	population = initialPopulation

	@assert popNum > 0
	@assert nMut >= 0
	@assert 0 <= mutationProb <= 1

	generationNumber = 0
	while !GAShouldStop(population, generationNumber)
		# Step 2
		fitness = map(GAFitness, population)
		constraintValues = LinearAlgebra.normalize!(map(
			entity -> map(constraint -> constraint(entity), constraints),
			population))

		# Equation 13
		constraintViolationValues = map(
			valueList -> sum(x -> max(0, x), valueList),
			constraintValues)

		# Equation 14
		numberOfViolations = map(
			valueList -> sum(x -> x > 0, valueList) / length(constraints),
			constraintValues)

		function customLessThan(entity1::Entity, entity2::Entity)::Bool
			entity1Index = findfirst(isequal(entity1), population)
			entity2Index = findfirst(isequal(entity2), population)
			entity1Feasible = numberOfViolations[entity1Index] == 0
			entity2Feasible = numberOfViolations[entity2Index] == 0

			if entity1Feasible && entity2Feasible
				# Winner is the one with the highest fitness value
				return fitness[entity1Index] > fitness[entity2Index]
			end

			if xor(entity1Feasible, entity2Feasible)
				# Winner is the feasible one
				return entity1Feasible
			end

			# Both are infeasible
			if numberOfViolations[entity1Index] == numberOfViolations[entity2Index]
				# Winner is the one with the lowest CV
				return constraintViolationValues[entity1Index] < constraintViolationValues[entity2Index]
			else
				# Winner is the one with the lowest NV
				return numberOfViolations[entity1Index] < numberOfViolations[entity2Index]
			end
		end

		# Step 3
		sortedPopulation = sort(population, lt=customLessThan)

		# Step 4
		elites = sortedPopulation[1:nElite]

		# Step 5
		crossed = map(_ -> GACrossover(rand(elites), rand(elites)), 1:nCross)
		mutated = map(_ -> GAMutate(rand(elites), mutationProb), 1:nMut)
		population = vcat(elites, crossed, mutated)
		generationNumber += 1
	end

	return population
end

global (limb, motors, gearRatios) = loadProblem(
	"res/constraints2.json",
	"HephaestusArmLimbOne",
	"res/motorOptions.json")

function GAFitness(entity::Entity)::Float64
	# Use negative of price because we maximize fitness
	return -1 * (motors[entity.motor1Index].price + motors[entity.motor2Index].price +
		motors[entity.motor3Index].price)
end

function GACrossover(entity1::Entity, entity2::Entity)::Entity
	# Take the motors  and ratios from entity1 and the links from entity2
	return Entity(entity1.motor1Index, entity1.motor2Index, entity1.motor3Index,
		entity2.link1Length, entity2.link2Length, entity2.link3Length,
		entity1.gearRatio1, entity1.gearRatio2, entity1.gearRatio3)
end

function GAMutate(entity::Entity, mutationProb::Float64)::Entity
	return Entity(
		if (rand() <= mutationProb) rand(1:length(motors)) else entity.motor1Index end,
		if (rand() <= mutationProb) rand(1:length(motors)) else entity.motor2Index end,
		if (rand() <= mutationProb) rand(1:length(motors)) else entity.motor3Index end,
		if (rand() <= mutationProb) rand(limb.minLinks[1].dhParam.r:limb.maxLinks[1].dhParam.r)/1000 else entity.link1Length end,
		if (rand() <= mutationProb) rand(limb.minLinks[2].dhParam.r:limb.maxLinks[2].dhParam.r)/1000 else entity.link2Length end,
		if (rand() <= mutationProb) rand(limb.minLinks[3].dhParam.r:limb.maxLinks[3].dhParam.r)/1000 else entity.link3Length end,
		if (rand() <= mutationProb) rand(gearRatios) else entity.gearRatio1 end,
		if (rand() <= mutationProb) rand(gearRatios) else entity.gearRatio2 end,
		if (rand() <= mutationProb) rand(gearRatios) else entity.gearRatio3 end
	)
end

function GAShouldStop(population::Vector{Entity}, generationNumber::Int64)::Bool
	return generationNumber > 200000
end

function makeRandomEntity()
	return GAMutate(Entity(0, 0, 0, 0, 0, 0, 0, 0, 0), 1.0)
	# mtrIndex = findfirst(x -> x.name == "stepperMotor-GenericNEMA14", motors)
	# return Entity(mtrIndex, mtrIndex, mtrIndex, 150, 150, 150, 0.05263157894736842, 0.05263157894736842, 0.05263157894736842)
end

function makeConstraints()
	return [entity::Entity -> (limb.tipForce * (entity.link1Length + entity.link2Length + entity.link3Length) +
				gravity * (motors[entity.motor2Index].mass * entity.link1Length +
				motors[entity.motor3Index].mass * (entity.link1Length + entity.link2Length))) -
				motors[entity.motor1Index].τStall / entity.gearRatio1,

		entity::Entity -> (limb.tipForce * (entity.link2Length + entity.link3Length) +
				motors[entity.motor3Index].mass * gravity * entity.link2Length) -
				motors[entity.motor2Index].τStall / entity.gearRatio2,

		entity::Entity -> (limb.tipForce * entity.link3Length) -
				motors[entity.motor3Index].τStall / entity.gearRatio3,

		entity::Entity -> (limb.tipVelocity / (entity.link1Length + entity.link2Length + entity.link3Length)) -
				motors[entity.motor1Index].ωFree * entity.gearRatio1,

		entity::Entity -> (limb.tipVelocity / (entity.link2Length + entity.link3Length)) -
				motors[entity.motor2Index].ωFree * entity.gearRatio2,

		entity::Entity -> (limb.tipVelocity / entity.link3Length) -
				motors[entity.motor3Index].ωFree * entity.gearRatio3,

		entity::Entity -> abs(entity.link1Length + entity.link2Length + entity.link3Length - 0.4) - 1e-4]
end

global finalPopulation = geneticAlgorithm(
	map(x -> makeRandomEntity(), 1:100),
	makeConstraints(),
	0.05,
	1,
	94)

global bestEntity = finalPopulation[argmax(map(GAFitness, finalPopulation))]
println("Motor 1=", motors[bestEntity.motor1Index],
	"\nMotor 2=", motors[bestEntity.motor2Index],
	"\nMotor 3=", motors[bestEntity.motor3Index],
	"\nLink 1=", bestEntity.link1Length * 1000,
	"\nLink 2=", bestEntity.link2Length * 1000,
	"\nLink 3=", bestEntity.link3Length * 1000,
	"\nGear ratio 1=", bestEntity.gearRatio1,
	"\nGear ratio 2=", bestEntity.gearRatio2,
	"\nGear ratio 3=", bestEntity.gearRatio3)

println("Fitness=", GAFitness(bestEntity))

global bestConstraintValues = map(x -> x(bestEntity), makeConstraints())
println("Constraint values: ", bestConstraintValues)
if max(bestConstraintValues...) <= 0
	println("Feasible.")
else
	println("Not feasible.")
end
