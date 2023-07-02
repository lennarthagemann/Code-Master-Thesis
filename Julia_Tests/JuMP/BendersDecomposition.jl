using JuMP, CPLEX
#=
 In Benders Decomposition we iteratively solve a mixed-integer linear program by solving smaller subproblems.
 We divide a problem by sorting out integer and continous variables, then solving the parametrized easier linear program.
 The linear program is convex with respect to variation in the integral part. Thus, we can construct subgradients for the 
 parametrized problem to add to the formulation. We reiterate until or close to optimality.
=#

# Function Parameters for Optimization

MAXIMUM_ITERATIONS = 100
ABSOLUTE_OPTIMALITY_GAP = 1e-6

#Simple Benders Problem data
c_1 = [1, 4]
c_2 = [2, 3]
dim_x = length(c_1)
dim_y = length(c_2)
b = [-2; -3]
A_1 = [1 -3; -1 -3]
A_2 = [1 -2; -1 -1]
M = -1000;
#Initialize model
model = Model(CPLEX.Optimizer)
set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-8)
@variable(model, x[1:dim_x] >= 0, Int)
@variable(model, θ >= M)
@objective(model, Min, c_1' * x + θ)
print(model)

function solve_subproblem(x)
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    @variable(model, y[1:dim_y] >= 0)
    con = @constraint(model, A_2 * y .<= b - A_1 * x)
    @objective(model, Min, c_2' * y)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (obj = objective_value(model), y = value.(y), π = dual.(con))
end

function print_iteration(k, args...)
    f(x) = JuMP.Printf.@sprintf("%12.4e", x)
    println(lpad(k, 9), " ", join(f.(args), " "))
    return
end

println("Iteration  Lower Bound  Upper Bound          Gap")
for k in 1:MAXIMUM_ITERATIONS
    optimize!(model)
    lower_bound = objective_value(model)
    x_k = value.(x)
    ret = solve_subproblem(x_k)
    upper_bound = c_1' * x_k + ret.obj
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution")
        break
    end
    cut = @constraint(model, θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    @info "Adding the cut $(cut)"
end

optimize!(model)
x_optimal = value.(x)

optimal_ret = solve_subproblem(x_optimal)
y_optimal = optimal_ret.y