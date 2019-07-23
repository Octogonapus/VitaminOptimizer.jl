"""
	optimalIndices(slots)

Find the indices of the chosen values of `slots` (the indices where `slots[i] ≈ 1`).
"""
optimalIndices(slots) = [findfirst(x -> abs(x - 1) <= 1e-5, value.(slot)) for slot in slots]

"""
	optimalColumns(featureMatrix::FeatureMatrix)

Get the columns of the `featureMatrix` corresponding to the chosen values for the slots.
"""
function optimalColumns(featureMatrix::FeatureMatrix)
	indices = optimalIndices(featureMatrix.slots)

	# Check if any indices are nothing because that will throw an exception later on
	for (index, slot) in zip(indices, featureMatrix.slots)
		if index == nothing
			firstNonzero = findall(x -> x != 0, value.(slot))

			for x in firstNonzero
				println(string(value.(slot)[x]))
			end

			throw(ErrorException("`nothing` in indices. First nonzero index: " * string(firstNonzero)))
		end
	end

	return [featureMatrix.matrix[:, i] for i in indices]
end

"""
	findMotorIndex(featureMatrixColumn, motors)

Find the index of the motor in the `motors` array by searching for a motor with τStall, ωFree, price,
and mass equal to those in the `featureMatrixColumn` (after un-applying the gear ratio).
"""
function findMotorIndex(featureMatrixColumn, motors)
	ratio = featureMatrixColumn[5]
	index = findfirst(
		# Approximate equality on τStall and ωFree because we are un-applying the gear ratio.
		x::Motor -> x.τStall ≈ featureMatrixColumn[1] * ratio &&
			x.ωFree ≈ featureMatrixColumn[2] / ratio &&
			x.price == featureMatrixColumn[3] &&
			x.mass == featureMatrixColumn[4],
		motors)

	if index == nothing
		throw(ErrorException("`nothing` as motor index. motors:\n" * string(motors)))
	end

	return index
end
