
"""	failedToOptimize(model)

Check if the `model` failed to optimize.
"""
failedToOptimize(model) = !(termination_status(model) == MOI.OPTIMAL ||
	(termination_status(model) == MOI.TIME_LIMIT && has_values(model)))
