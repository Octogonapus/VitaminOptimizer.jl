include("parseConstraints.jl")
include("parseMotorOptions.jl")
include("problemUtil.jl")

import LinearAlgebra
using Plots, Statistics, Distributed

const gravity = 9.80665

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

function Base.show(io::IO, entity::Entity)
	println(io,
		"Motor 1=", motors[entity.motor1Index],
		"\nMotor 2=", motors[entity.motor2Index],
		"\nMotor 3=", motors[entity.motor3Index],
		"\nLink 1=", entity.link1Length * 1000,
		"\nLink 2=", entity.link2Length * 1000,
		"\nLink 3=", entity.link3Length * 1000,
		"\nGear ratio 1=", entity.gearRatio1,
		"\nGear ratio 2=", entity.gearRatio2,
		"\nGear ratio 3=", entity.gearRatio3,
		"\nFitness=", GAFitness(entity),
		"\nConstraint values=", map(x -> x(entity), makeConstraints()),
		"\nFeasible=", isFeasible(entity))
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
	mutationProb::Float64, eliteRatio::Float64, crossRatio::Float64)
	popNum::Int64 = length(initialPopulation)
	nElite::Int64 = round(popNum * eliteRatio)
	nCross::Int64 = round(popNum * crossRatio)
	nMut::Int64 = popNum - nElite - nCross
	population = initialPopulation

	@assert popNum > 0
	@assert nMut >= 0
	@assert 0 <= mutationProb <= 1

	avgFitness::Array{Float64, 1} = []

	generationNumber = 0
	while !GAShouldStop(population, generationNumber)
		# Step 2
		fitness = map(GAFitness, population)
		push!(avgFitness, mean(fitness))

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

		function entityLessThan(entity1::Entity, entity2::Entity)::Bool
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
		sortedPopulation = sort(population, lt=entityLessThan)

		# Step 4
		elites = sortedPopulation[1:nElite]

		# Step 5
		crossed = map(_ -> GACrossover(rand(elites), rand(elites)), 1:nCross)
		mutated = map(_ -> GAMutate(rand(elites), mutationProb), 1:nMut)
		population = vcat(elites, crossed, mutated)
		generationNumber += 1

		if generationNumber % 1000 == 0
			println("Generation ", generationNumber)
		end
	end

	return population, avgFitness
end

global (limb, motors, gearRatios) = loadProblem(
	"res/constraints2.json",
	"HephaestusArmLimbOne",
	"res/motorOptions.json",
	1.0)

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
	return generationNumber > maxNumGenerations
end

"""
Each of these constraints is of the form g(x) <= 0.
"""
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

		entity::Entity -> abs(entity.link1Length + entity.link2Length + entity.link3Length - 0.4) - 1e-4,

		entity::Entity -> makeNLForceConstaint(entity)]
end

function makeNLForceConstaint(entity::Entity)::Float64
	θ3 = π/2 - atan(limb.targetZ/1000, limb.targetX/1000)
	h = sqrt((limb.targetX/1000)^2 + ((limb.targetZ/1000) - entity.link1Length)^2)

	if h > entity.link2Length + entity.link3Length
		return 1.0
	end

	θ1asinArg = (entity.link1Length * sin(θ3)) / h
	if abs(θ1asinArg) > 1
		return 1.0
	end
	θ1 = asin(θ1asinArg)

	βacosArg = (-h^2 + entity.link2Length^2 + entity.link3Length^2) / (2 * entity.link2Length * entity.link3Length)
	if abs(βacosArg) > 1
		return 1.0
	end
	β = acos(βacosArg)

	θ4asinArg = (entity.link2Length * sin(β)) / h
	if abs(θ4asinArg) > 1
		return 1.0
	end
	θ4 = asin(θ4asinArg)

	α = (π - β - θ4) + (π - θ1 - θ3)
	J = [0 (entity.link3Length * cos(α + β) - entity.link2Length * sin(α)) (entity.link3Length * cos(α + β));
	     (entity.link1Length + entity.link3Length * sin(α + β) + entity.link2Length * cos(α)) 0 0;
		 0 (-entity.link3Length * sin(α + β) - entity.link2Length * cos(α)) (-entity.link3Length * sin(α + β))]
	jointTorques = transpose(J) * [0, 0, limb.tipForceClosePos]

	# println(limb.tipForce)
	# println(entity.link1Length)
	# println(entity.link2Length)
	# println(entity.link3Length)
	# println(α)
	# println(β)
	# println(jointTorques)

	return (jointTorques[2] - motors[entity.motor2Index].τStall / entity.gearRatio2) +
	       (jointTorques[3] - motors[entity.motor3Index].τStall / entity.gearRatio3)
end

function makeRandomEntity()
	return GAMutate(Entity(0, 0, 0, 0, 0, 0, 0, 0, 0), 1.0)
end

function isFeasible(entity::Entity)
	constraintValues = map(x -> x(entity), makeConstraints())
	return maximum(constraintValues) <= 0
end

const global maxNumGenerations = 150000

global (finalPopulation, avgFitness) = geneticAlgorithm(
	map(x -> makeRandomEntity(), 1:100),
	makeConstraints(),
	0.05,
	0.25,
	1 - 0.25 - 0.05)

println("Final population:")
for x::Entity in finalPopulation
	println(x)
end

global feasibleEntities = filter(x -> isFeasible(x), finalPopulation)

if !isempty(feasibleEntities)
	global bestFitness = maximum(map(x -> GAFitness(x), feasibleEntities))
	global bestEntities = Set(filter(x -> GAFitness(x) ≈ bestFitness, feasibleEntities))

	println("Best entities:")
	for x::Entity in bestEntities
		println(x)
	end
end

avgFitnessPerGen = plot(1:(maxNumGenerations+1), avgFitness, title="Average Fitness per Generation",
	label=["Average Fitness"], xlabel="Generation Number")
savefig(avgFitnessPerGen, "average_fitness_per_generation.png")
