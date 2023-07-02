using JuMP, BilevelJuMP, CPLEX

model = BilevelModel(CPLEX.Optimizer, mode = BilevelJuMP.IndicatorMode())

@variable(Lower(model), x)
@variable(Upper(model), y, Int)

@objective(Upper(model), Min, 3x + y)
@constraints(Upper(model), begin
    x <= 5
    y <= 8
    y >= 0
end)

@objective(Lower(model), Min, -x)
@constraints(Lower(model), begin
     x +  y <= 8
    4x +  y >= 8
    2x +  y <= 13
    2x - 7y <= 0
end)

optimize!(model)

objective_value(model) # = 3 * (3.5 * 8/15) + 8/15 # = 6.13...
value(x) # = 3.5 * 8/15 # = 1.86...
value(y) # = 8/15 # = 0.53...