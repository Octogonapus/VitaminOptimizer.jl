"""
	loadProblem(constraintsFile::String, limbName::String, motorOptionsFile::String)

Load the constraints from `constraintsFile` and the motor options from
`motorOptionsFile`. Select limb `limbName` from the constraints.
"""
function loadProblem(constraintsFile::String, limbName::String, motorOptionsFile::String, linearConversion::Float64)
	limb = parseConstraints!(constraintsFile, [limbName], linearConversion)[1]
	motors = parseMotorOptions!(motorOptionsFile)

	ratios = collect(range(1, stop=limb.maxGearRatio, step=limb.gearRatioStep))
	gearRatios = Set(hcat(ratios, 1 ./ ratios))

	return (limb, motors, gearRatios)
end
