using Gurobi

function makeGurobiModel!(maxNumSolutions::Int64)::Model
	env = Gurobi.Env()
	setparam!(env, "PoolSearchMode", 2)
	setparam!(env, "PoolSolutions", 10)
	return Model(with_optimizer(Gurobi.Optimizer, env, Presolve=1))
end
