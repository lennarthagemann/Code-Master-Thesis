"""
Test File to create new features for the water regulation process or the optimizaiton models
"""

using JuMP
using CPLEX
using Distributions
using LinearAlgebra
using Statistics
using Dates
using DataFrames
using SDDP
using Plots
import CSV
using JSON
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet(pwd() * "\\Water_Regulation\\WaterRegulation.jl")
using .WaterRegulation

filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
filepath_prices = pwd() *  "\\Data\\Spot Prices\\prices_df.csv"
filepath_results = pwd() * "\\Results\\LambdaZero\\"
res, plants, parts = read_data(filepath_Ljungan)


k = plants[1]
ek = k.equivalent
Qspill = k.spill_reference_level
n = 5

"""
Calculate cuts for generation function. Given a number of cuts, efficiency e and spill limit, return coefficients of the cuts 
by linear interpolation of the points.
"""

kx, keffs, x, y = Generation_Cuts(Qspill, ek, n)

abstract type AbstractConfiguration end

struct AnticipatoryConfig <: AbstractConfiguration end

struct NonAnticipatoryConfig <: AbstractConfiguration end

function add_state_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, res_real_initial::Dict{Reservoir, Float64},res_ind_initial::Dict{Reservoir, Float64}, ::AnticipatoryConfig)
    @variables(subproblem, begin
        0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
        res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
        Qnom[r = res], (SDDP.State, initial_value = 0)
        d[t = 1:T], (SDDP.State, initial_value = 0)
    end)
    return
end

function add_state_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, res_real_initial::Dict{Reservoir, Float64}, res_ind_initial::Dict{Reservoir, Float64}, ::NonAnticipatoryConfig)
    @variables(subproblem, begin
        0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
        res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
        Qnom[r = res], (SDDP.State, initial_value = 0)
        d[t = 1:T], (SDDP.State, initial_value = 0)
    end)
    return
end

function add_control_variables(subproblem::Model, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, plants_O::Vector{HydropowerPlant}, T::Int64, ::AnticipatoryConfig)
    @variables(subproblem, begin
        Qnom_change[r = res] >= 0
        d_bid[t = 1:T] >= 0
        Qeff[k = plants_j, t = 1:T] >= 0
        Qreal[r = res, t = 1:T] >= 0
        Qadj[r = res] >= 0
        P_Swap[r = res]
        P_Over[k = plants_O] >= 0
        BALANCE_INDICATOR[r = res], Bin
        z_up[t = 1:T] >= 0
        z_down[t = 1:T] >= 0
    end)
    return
end

function add_control_variables(subproblem::Model, res::Vector{Reservoir}, plants_j:: Vector{HydropowerPlant}, T::Int64, ::NonAnticipatoryConfig)
    @variables(subproblem, begin
        d_bid[t = 1:T] >= 0
        Qeff[k = plants_j, t = 1:T] >= 0
        Qreal[r = res, t = 1:T] >= 0
        Qnom_change[r = res] >= 0
        BALANCE_INDICATOR[r = res], Bin
        z_up[t = 1:T] >= 0
        z_down[t = 1:T] >= 0
    end)
    return
end

function add_random_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, ::NonAnticipatoryConfig)
    @variables(subproblem, begin
        c[t = 1:T]
        Qinflow[r = res] >= 0
    end)
    return
end

function add_random_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, ::AnticipatoryConfig)   
    @variables(subproblem, begin
        c[t = 1:T]
        Qinflow[r = res] >= 0
        Qnom_O[r = res]
    end)
    return
end

function add_stage_objective(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, j::Participant, mean_price::Dict{Reservoir, Float64}, T::Int64,::NonAnticipatoryConfig)
    d = subproblem[:d]
    c = subproblem[:c]
    res_ind = subproblem[:res_ind]
    if node == 1
        @stageobjective(subproblem, 0)
    else
        @stageobjective(subproblem,  -(sum(c[t] * d[t].in for t in 1:T for k in plants_j)
        # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
        + sum((res_ind[r].out - res_ind[r].in) * j.participationrate[r] * mean_price[r] for r in res))/1e3)
    end
    return
end

function add_stage_objective(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, j::Participant, mean_price::Dict{Reservoir, Float64}, T::Int64, ::AnticipatoryConfig)
    d = subproblem[:d]
    c = subproblem[:c]
    res_ind = subproblem[:res_ind]
    P_Swap = subproblem[:P_Swap]
    if node == 1
        # Objective function
        @stageobjective(subproblem, 0)
    else
        # Objective Function
        @stageobjective(subproblem,  -(sum(c[t] * d[t].in for t in 1:T for k in plants_j)
        # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
        + sum((res_ind[r].out - res_ind[r].in)  * j.participationrate[r] * mean_price[r] for r in res))/1e3)
    end
    return
end

function add_transition_function(subproblem::Model, res::Vector{Reservoir}, node::Int64, T::Int64, Qref::Dict{Reservoir, Float64}, ::NonAnticipatoryConfig)
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    Qnom_change = subproblem[:Qnom_change]
    Qinflow = subproblem[:Qinflow]
    d = subproblem[:d]
    d_bid = subproblem[:d_bid]
    if node == 1
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in)
            @constraint(subproblem, res_ind[r].out == res_ind[r].in)
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
            @constratin(subproblem, d[t].out == d_bid[t])
        end
    else
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
            @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
            @constratin(subproblem, d[t].out == d_bid[t])
        end
    end
    return
end

function add_transition_function(subproblem::Model, res::Vector{Reservoir}, node::Int64, T::Int64, Qref::Dict{Reservoir, Float64}, ::AnticipatoryConfig)
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    Qnom_change = subproblem[:Qnom_change]
    Qadj = subproblem[:Qadj]
    Qinflow = subproblem[:Qinflow]
    if node == 1
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in)
            @constraint(subproblem, res_ind[r].out == res_ind[r].in)
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
            @constraint(subproblem, d[t].out == d_bid[t])
        end
    else
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qadj[r] - Qinflow[r]))
            @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
            @constraint(subproblem, d[t].out == d_bid[t])
        end
    end
end

function add_stage_constraints(subproblem::Model, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, Qref::Dict{Reservoir, Float64}, T::Int64, stage_count::Int64, node::Int64, ::NonAnticipatoryConfig)
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
    Qnom_change = subproblem[:Qnom_change]
    Qeff = subproblem[:Qeff]
    Qreal = subproblem[:Qreal]
    d = subproblem[:d]
    z_up = subproblem[:z_up]
    z_down = subproblem[:z_down]
    if node == 1
        for r in res
            # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
            @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
            # Constraints
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
        end
    else
        for r in res
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qnom[r].in)
            # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
            @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
        end
        for k in plants_j
            # @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            @constraint(subproblem, BALANCE_INDICATOR[k.reservoir] => {sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level})
        end
        for t in 1:T
            for k in plants_j
                @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
            @constraint(subproblem, d[t].in == sum(Qeff[k,t] for k in plants_j) + z_up[t] - z_down[t])
        end
    end
    return
end

function add_stage_constraints(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, plants_O::Vector{HydropowerPlant}, j::Participant, O::Participant, T::Int64, Qref::Dict{Reservoir, Float64}, stage_count::Int64, ::AnticipatoryConfig)
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
    Qnom_change = subproblem[:Qnom_change]
    Qeff = subproblem[:Qeff]
    Qreal = subproblem[:Qreal]
    Qadj = subproblem[:Qadj]
    P_Over = subproblem[:P_Over]
    P_Swap = subproblem[:P_Swap]
    Qnom_O = subproblem[:Qnom_O]
    d = subproblem[:d]
    z_up = subproblem[:z_up]
    z_down = subproblem[:z_down]
    if node == 1
        for k in plants_j
            # @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            @constraint(subproblem, BALANCE_INDICATOR[k.reservoir] => {sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level})
        end
        for r in res
            # Constraints
            # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
            @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
        end
    else
        for r in res
            # Constraints
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            @constraint(subproblem, Qadj[r] == (Qnom_O[r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
            @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
            # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
            @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
        end
        for k in plants_O
            @constraint(subproblem, P_Over[k] >=  (sum(Qadj[r] for r in find_us_reservoir(k.reservoir)) - k.spill_reference_level) * k.equivalent)
            # @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            @constraint(subproblem, BALANCE_INDICATOR[k.reservoir] => {sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level})
        end
        for t in 1:T
            for k in plants_j
                @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
            @constraint(subproblem, d[t] == sum(Qeff[k,t] * k.equivalent for k in plants_j) + sum(P_Swap[r] for r in res) + z_up - z_down)
        end
    end
    return
end

function fix_random_variables(subproblem, res, node, T, Ω, P, ::NonAnticipatoryConfig)
    c = subproblem[:c]
    Qinflow = subproblem[:Qinflow]
    if !(node == 1)
        SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
            for t in 1:T
                JuMP.fix(c[t], ω.price[t])
            end
            for r in res
                JuMP.fix(Qinflow[r], ω.inflow; force=true)
            end
        end
    end
    return
end

function fix_random_variables(subproblem, res, node, T, Ω, P, ::AnticipatoryConfig)
    c = subproblem[:c]
    Qinflow = subproblem[:Qinflow]
    Qnom_O = subproblem[:Qnom_O]
    if !(node == 1)
        SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
            for t in 1:T
                JuMP.fix(c[t], ω.price[t])
            end
            for r in res
                JuMP.fix(Qinflow[r], ω.inflow; force=true)
                try
                    JuMP.fix(Qnom_O[r], ω.nomination[r], force=true)
                catch
                    JuMP.fix(Qnom_O[r], ω.nomination[node][r], force=true)
                end
            end
        end 
    end
    return
end

function ShortTermOptimizationNoAnticipationDevelopment(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    scenario_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Ω,
    P,
    Qref,
    mean_price,
    price_sample;
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    optimizer = CPLEX.Optimizer,
    BIG_M = 5e4,
    stall_bound = SDDP.BoundStalling(5, 1e-2),
    config = NonAnticipatoryConfig()
    )
    res = filter(r -> j.participationrate[r] > 0, all_res)
    function subproblem_builder(subproblem::Model, node::Int)
        # Define State, Control and Random Variables
        add_state_variables(subproblem, res, T, res_real_initial, res_ind_initial, config)
        add_control_variables(subproblem, res, plants_j, T, config)
        add_random_variables(subproblem, res, T, config)
        # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
        add_transition_function(subproblem, res, node, T, Qref, config)
        # Add constraints for every stage
        add_stage_constraints(subproblem, res, plants_j, Qref, T, stage_count, node, config)
        # Fix Random Variables for nondeterministic stages
        fix_random_variables(subproblem, res, node, T, Ω, P, config)
        # Add Objective Function 
        add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
        return subproblem
    end
    # Define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
        optimizer = optimizer
    )
    # Train the model
    SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [stall_bound], print_level = printlevel, risk_measure = riskmeasure)
    # obtain decision rule in all steps
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:scenario_count)], inflow = 0),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        for r in filter(x -> !(x in res), all_res)
            Qnom[r] = Qref[r]
        end
        append!(rules, [rule])
        append!(nominations, [Qnom])
    end
    return model, rules, nominations
end

function ShortTermOptimizationAnticipationDevelopment(
    res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_O::Array{HydropowerPlant},
    parts::Array{Participant},
    plants_j::Array{HydropowerPlant},
    ITERATION_COUNT::Int64,
    stage_count::Int64,
    SCENARIO_COUNT::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Ω,
    P,
    nom,
    Qref,
    mean_price,
    price_sample;
    riskmeasure = SDDP.Expectation(),
    optimizer = CPLEX.Optimizer,
    printlevel = 1,
    BIG_M = 5e4,
    stall_bound = SDDP.BoundStalling(5, 1e-2),
    config = AnticipatoryConfig()
    )
    function subproblem_builder(subproblem::Model, node::Int)
        # Add State, Control and Random Variables
        add_state_variables(subproblem, res, T, res_real_initial, res_ind_initial, config)
        add_control_variables(subproblem, res, plants_j, plants_O, T, config)
        add_random_variables(subproblem, res, T, config)
        # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
        add_transition_function(subproblem, res, node, T, Qref, config)
        # Add constraints for every stage
        add_stage_constraints(subproblem, node, res, plants_j, plants_O, j, O, Qref, T, stage_count, config)
        # Fix Random Variables for nondeterministic stages
        fix_random_variables(subproblem, res, node, T, Ω, P, config)
        # Add Objective Function 
        add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
        return subproblem
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
        optimizer = optimizer
    )
    # train the model
    SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound], risk_measure = riskmeasure, print_level = printlevel)
    # obtain decision rule in all steps, as well as nominations for the reservoir situation.
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 0.1, nomination = nom[rand(1:SCENARIO_COUNT)]),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        push!(rules, rule)
        push!(nominations, Qnom)
    end
    return model, rules, nominations
end

# model, rules, nominations = ShortTermOptimizationNoAnticipationDevelopment(
#     res,
#     parts[1],
#     OtherParticipant(parts[1], parts),
#     parts[1].plants,
#     10,
#     7,
#     3,
#     24,

# )

function MediumTermOptimizationDevelopment(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    scenario_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Ω,
    P,
    Qref,
    mean_price,
    price_sample;
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    optimizer = CPLEX.Optimizer,
    BIG_M = 5e4,
    stall_bound = SDDP.BoundStalling(5, 1e-2),
    config = NonAnticipatoryConfig()
)

    return
end