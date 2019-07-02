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
motorNames = keys(motorData)
F_m = hcat([
	[motorData[mtr]["MaxTorqueNewtonmeters"],
	 motorData[mtr]["MaxFreeSpeedRadPerSec"],
	 motorData[mtr]["price"],
	 motorData[mtr]["massKg"]]
	for mtr in motorNames]...)

# Select first two cols to keep it simple for now
# F_m = F_m[:, 1:2]

gravity = 9.81

env = Gurobi.Env()
setparam!(env, "LogFile",
		  "/home/salmon/Documents/auto-configured-vitamins-optimizer/main.log")

# pass params as keyword arguments to GurobiSolver
model = Model(with_optimizer(Gurobi.Optimizer, Presolve=1))
# model = Model(with_optimizer(GLPK.Optimizer))

(numRows, numCols) = size(F_m)

τRow = [1 0 0 0]
ωRow = [0 1 0 0]
priceRow = [0 0 1 0]
massRow = [0 0 0 1]

# Each slot is a binary vector with a 1 that picks which motor to use.
@variable(model, slot1[1:numCols], Bin)
@constraint(model, slot1Unique, sum(slot1) == 1)

@variable(model, slot2[1:numCols], Bin)
@constraint(model, slot2Unique, sum(slot2) == 1)

@variable(model, slot3[1:numCols], Bin)
@constraint(model, slot3Unique, sum(slot3) == 1)

# Equation 3
@constraint(model, eq3, τRow * F_m * slot1 .>=
                        tipForce * (linkDhA[1] + linkDhA[2] + linkDhA[3]) +
                        gravity * (massRow * F_m * slot2 * linkDhA[1] +
                        massRow * F_m * slot3 * (linkDhA[1] + linkDhA[2])))

# Equation 4
@constraint(model, eq4, τRow * F_m * slot2 .>=
                        tipForce * (linkDhA[2] + linkDhA[3]) +
                        massRow * F_m * slot3 * gravity * linkDhA[2])

# Equation 5
@constraint(model, eq5, τRow * F_m * slot3 .>=
                        tipForce * linkDhA[3])

@objective(model, Min, sum(x -> priceRow * F_m * x, [slot1, slot2, slot3])[1])

optimize!(model)

println("Optimal objective: ", objective_value(model),
	". slot1 = ", value.(slot1), ", slot2 = ", value.(slot2),
	", slot3 = ", value.(slot3))
