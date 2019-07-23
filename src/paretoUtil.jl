function stayAtParetoFrontier(model::Model, objectiveFunction)
	@constraint(model, objective_value(model) - objectiveFunction <= 1e-5)
end

"""
	exploreParetoFrontier(model::Model, featureMatrix::FeatureMatrix, motors, filename::String)

Iteratively optimize the `model` to find all solutions at a given `optimalObjectiveValue` by adding a
constraint to disallow the most recent combination of slot values. Returns an array of all solutions.
"""
function exploreParetoFrontier(model::Model, featureMatrix::FeatureMatrix, motors, filename::String)
	# Disallow the current solution by disallowing the combination of the current slot1, slot2, and slot3
	# values.
	@constraint(
		model,
		sum(x -> x[1][x[2]], zip(featureMatrix.slots, optimalIndices(featureMatrix.slots)))
			<= length(featureMatrix.slots) - 1)

	# Optimize again to find a different solution.
	optimize!(model)

	# If we are no longer at the given optimum or if the model failed to opimize, stop.
	if failedToOptimize(model)
		return []
	else
		# Record the current solution.
		solution = findOptimalChoices(featureMatrix, motors)
		printOptimizationResult!(model, solution)
		saveOptimizationResult(model, solution, filename)

		# Keep finding more solutions.
		otherSolutions = exploreParetoFrontier(model, featureMatrix, motors, filename)

		# Add the current solution to the end of the other solutions.
		if otherSolutions == []
			return solution
		else
			return vcat(solution, otherSolutions)
		end
	end
end
