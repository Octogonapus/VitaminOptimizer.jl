import JSON, GLPK
using JuMP, Gurobi

j = JSON.parsefile("res/constraints1.json")
tipVelocity = j["requiredTipVelocityMeterPerSec"]
tipForce = j["requiredTipForceNewtons"]
limb = j["HephaestusArmLimbOne"]
limbConfig = limb["min"]
linkNames = sort(collect(keys(limbConfig)), by=x -> limbConfig[x]["index"])
# Convert millimeters to meters
linkDhA = [limbConfig[name]["dh-A"] / 1000 for name=linkNames]

j = JSON.parsefile("res/motorOptions.json")
motorData = j["data"]
motorNames = collect(keys(motorData))
gearRatios = [7, 5, 3, 1, 1/3, 1/5, 1/7]
# Take the coproduct of motors with gear ratios to generate all possible combinations.
motorNamesWithRatios = [(name, ratio) for name in motorNames, ratio in gearRatios]
F_m = hcat([
	[motorData[mtr]["MaxTorqueNewtonmeters"] / ratio,
	 motorData[mtr]["MaxFreeSpeedRadPerSec"] * ratio,
	 motorData[mtr]["price"],
	 motorData[mtr]["massKg"],
	 (motorData[mtr]["MaxFreeSpeedRadPerSec"] * ratio) / (motorData[mtr]["MaxTorqueNewtonmeters"] / ratio)]
	for (mtr, ratio) in motorNamesWithRatios]...)

# Select first two cols to keep it simple for now
# F_m = F_m[:, 1:2]

const gravity = 9.80665

env = Gurobi.Env()
setparam!(env, "LogFile",
		  "/home/salmon/Documents/auto-configured-vitamins-optimizer/main.log")

# pass params as keyword arguments to GurobiSolver
model = Model(with_optimizer(Gurobi.Optimizer, Presolve=1))
# model = Model(with_optimizer(GLPK.Optimizer))

(numRows, numCols) = size(F_m)

τRow = [1 0 0 0 0]
ωRow = [0 1 0 0 0]
priceRow = [0 0 1 0 0]
massRow = [0 0 0 1 0]
omegaFuncRow = [0 0 0 0 1]

# Each slot is a binary vector with a 1 that picks which motor to use.
@variable(model, slot1[1:numCols], Bin)
@constraint(model, slot1Unique, sum(slot1) == 1)

@variable(model, slot2[1:numCols], Bin)
@constraint(model, slot2Unique, sum(slot2) == 1)

@variable(model, slot3[1:numCols], Bin)
@constraint(model, slot3Unique, sum(slot3) == 1)

allSlots = [slot1, slot2, slot3]

"""
	@ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i::Int64)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Array{GenericAffExpr{Float64,VariableRef},1}, i::Int64) =
	((τRow * F_m * allSlots[i]) - τ) * (omegaFuncRow * F_m * allSlots[i])

"""
	@ω(τ::Float64, i::Int64)

Calculate the rotational speed given the applied torque `τ` for the motor in
slot `i`.
"""
ω(τ::Float64, i::Int64) =
	((τRow * F_m * allSlots[i]) .- τ) * (omegaFuncRow * F_m * allSlots[i])

# Equation 3
# linkDhA[1] is 0 which makes this boring
@expression(model, τ1Required, tipForce * (linkDhA[1] + linkDhA[2] + linkDhA[3]) +
							   gravity * (massRow * F_m * slot2 * linkDhA[1] +
							   massRow * F_m * slot3 * (linkDhA[1] + linkDhA[2])))
@constraint(model, eq3, τRow * F_m * slot1 .>= τ1Required)

# Equation 4
@expression(model, τ2Required, tipForce * (linkDhA[2] + linkDhA[3]) +
							   massRow * F_m * slot3 * gravity * linkDhA[2])
@constraint(model, eq4, τRow * F_m * slot2 .>= τ2Required)

# Equation 5
@expression(model, τ3Required, tipForce * linkDhA[3])
@constraint(model, eq5, τRow * F_m * slot3 .>= τ3Required)

# Equation 6
@expression(model, ω1Required, tipVelocity / (linkDhA[1] + linkDhA[2] + linkDhA[3]))
@constraint(model, eq6, ωRow * F_m * slot1 .>= ω1Required)

# Equation 7
@expression(model, ω2Required, tipVelocity / (linkDhA[2] + linkDhA[3]))
@constraint(model, eq7, ωRow * F_m * slot2 .>= ω2Required)

# Equation 8
@expression(model, ω3Required, tipVelocity / linkDhA[3])
@constraint(model, eq8, ωRow * F_m * slot3 .>= ω3Required)

@objective(model, Min, sum(x -> priceRow * F_m * x, allSlots)[1])

optimize!(model)

println("Optimal objective: ", objective_value(model),
	". slot1 = ", value.(slot1), ", slot2 = ", value.(slot2),
	", slot3 = ", value.(slot3))

optimalMotorIndices = [findfirst(isequal(1), value.(slot)) for slot in allSlots]
optimalMotorNames = [motorNamesWithRatios[i] for i in optimalMotorIndices]
println("Optimal motors: ", optimalMotorNames)
