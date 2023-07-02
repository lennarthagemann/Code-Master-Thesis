using JuMP, CPLEX


MAXIMUM_ITERATIONS = 100
ABSOLUTE_OPTIMALITY_GAP = 1e-6

c_1 = [1, 4]
c_2 = [2, 3]
dim_x = length(c_1)
dim_y = length(c_2)
b = [-2; -3]
A_1 = [1 -3; -1 -3]
A_2 = [1 -2; -1 -1]
M = -1000;

model = Model(CPLEX.Optimizer)
@variable(model, x[1:dim_x] >= 0, Int)
@variable(model, θ >= M)
@objective(model, Min, c_1' * x + θ)
print(model)

subproblem = Model(CPLEX.Optimizer)
@variable(subproblem, x_copy[1:dim_x])
@variable(subproblem, y[1:dim_y] >= 0)
@constraint(subproblem, A_1 * x_copy + A_2 * y .<= b)
@objective(subproblem, Min, c_2' * y)
print(subproblem)

function solve_subproblem(model, x)
    fix.(model[:x_copy], x) #value of x and copy are supposed to be the same, so first and second stage have same solution
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (
        obj = objective_value(model),
        y = value.(model[:y]),
        π = reduced_cost.(model[:x_copy]), #valid subgradient with respect to x
    )
end

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(subproblem, x_k)
    upper_bound = c_1' * x_k + ret.obj
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + ret.π' * (x .- x_k)) # no need to multiply ret.π' with -A_1, beacuse it is a valid subgradient with respect to x already. (saves costs)
    @info "Adding the cut $(cut)"
end

optimize!(model)
x_optimal = value.(x)

optimal_ret = solve_subproblem(subproblem, x_optimal)
y_optimal = optimal_ret.y