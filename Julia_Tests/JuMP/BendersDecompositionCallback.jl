using JuMP, CPLEX

ABSOLUTE_OPTIMALITY_GAP = 1e-8

#Add new constraints !while! solver is running (lazy constraint callbacks)

c_1 = [1, 4]
c_2 = [2, 3]
dim_x = length(c_1)
dim_y = length(c_2)
b = [-2; -3]
A_1 = [1 -3; -1 -3]
A_2 = [1 -2; -1 -1]
M = -1000;

lazy_model = Model(CPLEX.Optimizer)
@variable(lazy_model, x[1:dim_x] >= 0, Int)
@variable(lazy_model, θ >= M)
@objective(lazy_model, Min, θ)
print(lazy_model)


k = 0

function solve_subproblem(x)
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    set_silent(model)
    @variable(model, y[1:dim_y] >= 0)
    con = @constraint(model, A_2 * y .<= b - A_1 * x)
    @objective(model, Min, c_2' * y)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    return (obj = objective_value(model), y = value.(y), π = dual.(con))
end

function print_iteration(k, args...) # Efficiently prints arguments in table format
    f(x) = JuMP.Printf.@sprintf("%12.4e", x) #@sprintf -> Macro to print output formatted as string, from Printf standard library, e stands for scientific notation
    println(lpad(k, 9), " ", join(f.(args), " ")) #lpad adds whitespace before, join adds iterator together with the format given at the end 
    return
end



"""
    my_callback(cb_data)

A callback that implements Benders decomposition.
"""
function my_callback(cb_data)
    global k += 1
    x_k = callback_value.(cb_data, x) #primal solution of variable inside callback
    θ_k = callback_value(cb_data, θ)
    lower_bound = c_1' * x_k + θ_k
    ret = solve_subproblem(x_k)
    upper_bound = c_1' * x_k + c_2' * ret.y
    gap = (upper_bound - lower_bound) / upper_bound
    print_iteration(k, lower_bound, upper_bound, gap)
    if gap < ABSOLUTE_OPTIMALITY_GAP
        println("Terminating with the optimal solution.")
        return
    end
    cut = @build_constraint(θ >= ret.obj + -ret.π' * A_1 * (x .- x_k))
    MOI.submit(lazy_model, MOI.LazyConstraint(cb_data), cut)
    return
end

MOI.set(lazy_model, MOI.LazyConstraintCallback(), my_callback)
set_silent(lazy_model)
optimize!(lazy_model)
JuMP.Printf.@printf("Optimal solution is as following: %.2f", value.(x))