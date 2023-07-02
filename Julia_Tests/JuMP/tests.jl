using JuMP, CPLEX

model = Model(CPLEX.Optimizer)
set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-8)

@variable(model, 0 <= x <= 2.5, Int)
@variable(model, 0 <= y <= 2.5, Int)
@objective(model, Max, y)

optimize!(model)

println("x : $(value.(x))")
println("y : $(value.(y))")
