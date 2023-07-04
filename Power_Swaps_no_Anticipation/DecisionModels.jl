using JuMP
using CPLEX
using Distributions
using LinearAlgebra
using Statistics
using DataFrames
using SDDP
using Plots
import CSV
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation

T = 24
STAGE_COUNT = 7
ITERATION_COUNT = 50
SCENARIO_COUNT = 10
BIG_M = 10e+4

filepath_systemA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/SimpleReservoirSystem.json"
res, plants, parts = read_data(filepath_systemA)

INITIAL_RESERVOIR = Dict{Reservoir, Float64}(r => r.currentvolume for r in res)
INITIAL_INDIVIDUAL_RESERVOIR = Dict{Reservoir, Float64}(r => r.currentvolume for r in res)

filepath_prices = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Data/Spot Prices/prices_df.csv"
price_data = CSV.read(filepath_prices, DataFrame)
price_data = coalesce.(price_data, 46.79)
# Get Array of every row in price_data (excluding Data Column1)
price_scenarios = []
for row in eachrow(select(price_data, Not(:Column1)))
    push!(price_scenarios, collect(values(row)))
end
price_sample = Dict(s => price_scenarios[rand(1:length(price_scenarios), scenario_count)] for s in 1:stage_count)
max_hourly_price = max([max([mean(s) for s in sample]...) for sample in values(price_sample)]...)
mean_hourly_price = mean(mean(mean(mean(sample, dims=2)) for sample in values(price_sample)))
inflow_scenarios = [0.1]
Ω = Dict(s => [(price = c, inflow = Q) for c in price_sample[s] for Q in inflow_scenarios] for s in 1:stage_count)
P = Dict(s => [1/length(eachindex(Ω[s])) for i in eachindex(Ω)] for s in 1:stage_count)
Qref = Dict{Reservoir, Float64}(r => mean(inflow_scenarios) for r in res)
# Plot each of the price_samples. 

"""
Test Funktion, klappt bisher nur mit einem Reservoir, nicht innerhalb einer Kaskade.
"""
function ShortTermOptimizationNoAnticipation(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    printlevel::Bool;
    Ω = Ω,
    P = P,
    Qref = Qref,
    optimizer = CPLEX.Optimizer
    )
    res = filter(r -> j.participationrate[r] > 0, all_res)
    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_initial[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_initial[r])
        @variable(subproblem, Qnom[r = res], SDDP.State, initial_value = 0)
        # Control Variables
        @variable(subproblem, Qeff[k = plants_j, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, 0 <= Qnom_change[r = res] <= 120)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        # Random Variables
        @variable(subproblem, c[t = 1:T])
        @variable(subproblem, Qinflow[r = res] >= 0)
        if node == 1
            for r in res
                # Transition Function
                @constraint(subproblem, res_real[r].out == res_real[r].in)
                @constraint(subproblem, res_ind[r].out == res_ind[r].in)
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                # Constraints
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            end
            # Objective function
            @stageobjective(subproblem, 0)
        end
        if (node in 2:stage_count)
            SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
                for t in 1:T
                    JuMP.fix(c[t], ω.price[t])
                end
                for r in res
                    JuMP.fix(Qinflow[r], ω.inflow; force=true)
                end
            end 
                # Transition Function and Constraints
            for r in res
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= (res_real[r].out))
                @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qnom[r].in)
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            for k in plants_j
                @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            end
            for t in 1:T
                for k in plants_j
                    @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                    @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
                end
            end
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
            + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price
            + sum(( 1 - (j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_ind[r].out - res_ind[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price)/1e3)
            return subproblem
        end
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_hourly_price * sum(k.equivalent for k in plants_j) for r in res)/1e3,
        optimizer = optimizer
    )
    # train the model
    if printlevel == true
        SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [stall_bound])
    else
        SDDP.train(model; iteration_limit = iteration_count, print_level = 0, stopping_rules = [stall_bound])
    end
    # obtain decision rule in first step
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:scenario_count)], inflow = 20),
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
"""
Diese Funktion ist nur logisch für ein System mit einem einzelnen Reservoir!
(Ich hasse diese blöde Kaskade, es könnte alles so schön sein)
"""
function ShortTermOptimizationAnticipation(
    res::Array{Reservoir},
    parts::Array{Participant},
    j::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Qref = Qref,
    Ω = Ω,
    P = [1/length(eachindex(Ω)) for i in eachindex(Ω)],
    optimizer = CPLEX.Optimizer
    )
    O, plants_O = OtherParticipant(parts, j, res)
    _,_, Qnom_O = ShortTermOptimizationNoAnticipation(res, j, O, plants_O, iteration_count,stage_count, T, res_real_initial, res_ind_initial, true)
    @assert length(Qnom_O) == stage_count
    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_initial[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_initial[r])
        @variable(subproblem, Qnom[r = res], SDDP.State, initial_value = 0)
        # Control Variables
        @variable(subproblem, Qnom_change[r = res] >= 0)
        @variable(subproblem, Qeff[k = plants_j, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, Qadj[r = res])
        @variable(subproblem, Qmax[k = plants_O] >= 0)
        @variable(subproblem, P_Swap[r = res])
        @variable(subproblem, P_Over[k = plants_O] >= 0)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        # Random Variables
        @variable(subproblem, c[t = 1:T])
        @variable(subproblem, Qinflow[r = res] >= 0)
        if node == 1
            for r in res
                # Transition Function
                @constraint(subproblem, res_real[r].out == res_real[r].in)
                @constraint(subproblem, res_ind[r].out == res_ind[r].in)
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                # Constraints
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
                for k in plants_O
                    @constraint(subproblem, Qnom_change[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
                end
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            # Objective function
            @stageobjective(subproblem, 0)
        end
        if (node in 2:stage_count)
            SDDP.parameterize(subproblem, Ω, P) do ω
                for t in 1:T
                    JuMP.fix(c[t], ω.price[t])
                end
                for r in res
                    JuMP.fix(Qinflow[r], ω.inflow; force=true)
                end
            end 
            for r in res
                # Transition Function
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qadj[r]- Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                # Constraints
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
                @constraint(subproblem, Qadj[r] == (Qnom_O[node-1][r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
                @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
                @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
                for k in plants_O
                    @constraint(subproblem, Qnom_change[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
                end
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            for k in plants_O
                @constraint(subproblem, P_Over[k] >= (Qadj[k.reservoir] - k.spill_reference_level) * k.equivalent)
            end
            for t in 1:T
                for k in plants_j
                    @constraint(subproblem, Qeff[k, t] <= Qreal[k.reservoir, t])
                    @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
                end
            end
            # Objective Function
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
            + sum(c[t] * P_Swap[r] for t in 1:T for r in res)
            + sum((1 - (j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price
            + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_ind[r].out - res_ind[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price))
            return subproblem
        end
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_hourly_price * sum(k.equivalent for k in plants_j) for r in res),
        optimizer = optimizer
    )
    # train the model
    SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [SDDP.BoundStalling(2, 1e-4)])
    # obtain decision rule in all steps, as well as nominations for the reservoir situation.
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:scenario_count)], inflow = 20),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        push!(rules, rule)
        push!(nominations, Qnom)
    end
    return model, rules, nominations
end
"""
Include the deterministic equivalent of the problem from other users into every subproblem for every scenario.
Then we can have different nomination scenarios for the other producer instead of a fixed value.
"""
function ShortTermOptimizationWithAnticipationDistribution(
    res::Array{Reservoir},
    parts::Array{Participant},
    j::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    T::Int64,
    res_real_level::Dict{Reservoir, Float64},
    res_ind_level::Dict{Reservoir, Float64},
    Qref::Dict{Reservoir, Float64},
    Ω = [(price = c, inflow = Q) for c in price_sample for Q in inflow_scenarios],
    P = [1/length(eachindex(Ω)) for i in eachindex(Ω)],
    Optimizer = CPLEX.Optimizer
    )

    O = OtherParticipant(parts, j, res)[1]
    plants_O = OtherParticipant(parts, j, res)[2]
    # Create the Distribution of other nominations first

    Other_Nomination_Scenario = Dict((price, inflow, stage) => DeterministicNoAnticipation(
            res,
            O,
            plants_O,
            5,
            1,
            T,
            map(x -> x + 0.1* stage, price),
            Dict(r => float(inflow) for r in res),
            Dict(r => res_real_level[r] - stage * Qref[r] for r in res),
            Dict(r => res_ind_level[r] for r in res),
            Qref
            )
        for (price, inflow) in Ω for stage in 1:stage_count
    )

    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_level[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_level[r])
        # Control Variables
        @variable(subproblem, Qnom[r = res])
        @variable(subproblem, Qeff[k = plants_j, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, Qadj[r = res])
        @variable(subproblem, Qmax[k = plants_O] >= 0)
        @variable(subproblem, P_Swap[r = res])
        @variable(subproblem, P_Over[k = plants_O] >= 0)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        # Random Variables
        @variable(subproblem, c[t = 1:T])
        @variable(subproblem, Qinflow[r = res] >= 0)
        @variable(subproblem, Qnom_O[r = res] >= 0)
        SDDP.parameterize(subproblem, Ω, P) do ω
            for t in 1:T
                for r in res
                    if typeof(r.upstream_reservoir) == Nothing
                        JuMP.fix(Qinflow[r], ω.inflow; force=true)
                    else
                        JuMP.fix(Qinflow[r], 0; force=true)
                    end
                end
                try
                    JuMP.fix(c[t], ω.price[t] + 0.1 * float(node))
                catch err
                    println(ω)
                    println("Something has gone wrong. The price is... $(ω.price), the node is $(node), and the time is $(t).")
                end
            end
            for r in res
                JuMP.fix(Qnom_O[r], Other_Nomination_Scenario[(ω.price, ω.inflow, node)][r]; force=true)
            end
        end 
            # Transition Function and Constraints
        for r in res
            @constraint(subproblem, Qadj[r] == (Qnom[r] * j.participationrate[r] + Qnom_O[r] * O.participationrate[r])/(j.participationrate[r] + O.participationrate[r]))
            if typeof(r.upstream_reservoir) == Nothing
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qadj[r] - Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r] - Qref[r]))
            else
                @assert typeof(r.upstream_reservoir) == Vector{Reservoir} "type of upstream_reservoir is $(typeof(r.upstream_reservoir))"
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * Qadj[r])
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r] - Qref[r]) + T*(sum((Qnom[ur] - Qref[ur]) for ur in r.upstream_reservoir)))
            end
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
            for k in plants_O
                @constraint(subproblem, Qnom[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, P_Over[k] >= (Qadj[k.reservoir] - k.spill_reference_level) * k.equivalent)
            end
            @constraint(subproblem, P_Swap[r] ==  (Qnom[r] - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
            @constraint(subproblem, Qnom[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, res_real[r].out >= 0.8 * res_real[r].in)
            # Environmental Constraints, or constraints that enforce a lower bound at the end of the time horizon, as to not make the VF mad.
            if node == stage_count
                @constraint(subproblem, res_real[r].out >= 0.2 * res_real_level[r])
            end
        end
        for t in 1:T
            for k in plants_j
                @constraint(subproblem, Qeff[k, t] <= Qreal[k.reservoir, t])
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
        end
        # Stage Objective
        @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j) + sum(c[t] * P_Swap[r] for t in 1:T for r in res) + sum((res_real[r].out - res_real[r].in)  * sum(k.equivalent for k in filter(x -> x.reservoir in find_ds_reservoirs(r), plants_j)) for r in res) * mean_hourly_price))
        return subproblem
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(res_real_level[r] * max_hourly_price * sum(k.equivalent for k in plants_j) for r in res),
        optimizer = Optimizer
    )
    SDDP.train(
        model;
        iteration_limit = iteration_count,
        stopping_rules = [SDDP.BoundStalling(2, 1e-4)]
    )
    rules = []
    for stage in 1:stage_count
        push!(rules, SDDP.DecisionRule(model; node = stage))
    end
    return model, rules
end
"""
Same as the ShortTermOptimizationNoAnticipation, but with a deterministic price and inflow.
It is used in the scenario generation of ShortTermOptimizationWithAnticipation: After the uncertainty is revealed, we assume that
the other user also decided in that fashion.
"""
function DeterministicNoAnticipation(
    res::Array{Reservoir},
    j::Participant,
    plants::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    T::Int64,
    c::Array{Float64},
    Qinflow::Dict{Reservoir, Float64},
    res_real_level::Dict{Reservoir, Float64},
    res_ind_level::Dict{Reservoir, Float64},
    Qref = Dict{Reservoir, Float64}(r => 20 for r in res);
    WaterValue::Bool = true
    )
    @assert length(c) == T "The time horizon of the price scenario is not equal to T."
    # create function for subproblem_builder
    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_level[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_level[r])
        # Control Variables
        @variable(subproblem, Qnom[r = res] >= 0)
        @variable(subproblem, Qeff[k = plants, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        # Transition Function and Constraints
        for r in res
            if typeof(r.upstream_reservoir) == Nothing
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r] - Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r] - Qref[r]))
            else
                @assert typeof(r.upstream_reservoir) == Vector{Reservoir} "type of upstream_reservoir is $(typeof(r.upstream_reservoir))"
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * Qnom[r])
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r] - Qref[r]) + T*(sum((Qnom[ur] - Qref[ur]) for ur in r.upstream_reservoir)))
            end
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qnom[r])
            for k in plants
                @constraint(subproblem, Qnom[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
            end
            @constraint(subproblem, Qnom[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
        end
        for t in 1:T
            for k in plants
                @constraint(subproblem, Qeff[k, t] <= Qreal[k.reservoir, t])
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
        end

        # Stage Objective
        if WaterValue
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants) + sum((res_real[r].out - res_real[r].in)  * sum(k.equivalent for k in filter(x -> x.reservoir in find_ds_reservoirs(r), plants)) for r in res) * mean_hourly_price))
        else
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants)))
        end
        return subproblem
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * max_hourly_price * sum(k.equivalent for k in plants) for r in res),
        optimizer = CPLEX.Optimizer
    )
    # train the model
    SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [SDDP.BoundStalling(2, 1e-4)], print_level = 0)
    # obtain decision rule in first step
    rule = SDDP.DecisionRule(model; node = 1)
    solution = SDDP.evaluate(
        rule;
        incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_level[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_level[r] for r in res)),
        controls_to_record = [:Qnom]
    )
    Qnom_det = Dict(r => solution.controls[:Qnom][r] for r in res)

    return Qnom_det
end

function ObtainNomination(
    rules::Dict{Participant, SDDP.DecisionRule{Int64}},
    parts::Array{Participant},
    INITIAL_STATE::Dict{Reservoir, Float64},
    noise_vector
    )
    @assert length(keys(rules)) == length(parts)
    Qnoms = Dict{Participant, Dict{Reservoir, Float64}}(p => Dict{Reservoir, Float64}(r => 0.0 for r in res) for p in parts)
    # Evaluate each optimization problem by the obtained rule and the inital state.
    evals = Dict(p => SDDP.evaluate(
        rules[p];
        incoming_state = merge(Dict(Symbol("res_real[$(r)]") => INITIAL_STATE[r] for r in res), Dict(Symbol("res_ind[$(r)]") => INITIAL_STATE[r] for r in res)),
        noise = noise_vector,
        controls_to_record = [:Qnom, :res_real],
    ) for p in parts)
    # ObtainNomination Qnom for every stage from these Evaluations.
    for p in parts
        Qnoms[p] = Dict(r => evals[p].controls[:Qnom][r] for r in res)
    end
    return Qnoms
end

function DisplayDecisions(parts::Array{Participant},  simulations, variables, stage_count = STAGE_COUNT)
    # Give a pretty display of all the  decisions for every stage of trained policies for each producer.
    for s in 1:stage_count
        println("Stage $(s): ")
        for p in parts
            for v in variables
                println("\n $(v): ")
                for k in keys(simulations[p][1][s][v])
                    print(simulations[p][1][s][v][k], " - ")
                end
            end
        end
        println()
    end
end

function SimulatePolicies(models::Dict{Participant, SDDP.PolicyGraph{Int64}}, parts::Array{Participant}, tracked_variables::Array{Symbol})
    simulations = Dict{Participant, Any}(p => SDDP.simulate(
        # Trained model to simulate
        models[p],
        # Number of replications
        100,
        # A list of names to record
        tracked_variables,
    ) for p in parts)
    return simulations
end



# Create a dictionary of participants and their plants
part_plants = Dict{Participant, Array{HydropowerPlant}}(p => p.plants for p in parts)
j = parts[1]

# Get a decision rule for every power producer by feeding their plants into the ShortTermOptimizationNoAnticipation
# Functions return model, rules, nominations
# individual_solutions_anticipatory = Dict(p => ShortTermOptimizationAnticipation(
#     res,
#     parts,
#     p,
#     part_plants[p],
#     ITERATION_COUNT,
#     STAGE_COUNT,
#     T,
#     INITIAL_RESERVOIR,
#     INITIAL_INDIVIDUAL_RESERVOIR
#     )
# for p in parts)

# individual_solutions_nonanticipatory = Dict(p => ShortTermOptimizationNoAnticipation(
#     res,
#     OtherParticipant(parts, p, res)[1],
#     p,
#     part_plants[p],
#     ITERATION_COUNT,
#     STAGE_COUNT,
#     T,
#     INITIAL_RESERVOIR,
#     INITIAL_INDIVIDUAL_RESERVOIR,
#     true
#     )
# for p in parts)

# for p in parts
#     println("Nomination by $(p.name) using the anticipatory model: ", individual_solutions_anticipatory[p][3])
# #    println("Nomination by $(p.name) using the nonanticipatory model: ",individual_solutions_nonanticipatory[p][3])
# end

# decision_rules = Dict(p => individual_solutions_anticipatory[p][2] for p in parts)
# # Use the decision_rules to obtain a nomination under the currrent conditions.
# Qnoms = Dict{Participant, Dict{Reservoir, Float64}}(p => Dict{Reservoir, Float64}(r => 0.0 for r in res) for p in parts)

# for p in parts
#     Qnoms[p] = Dict(r => SDDP.evaluate(
#         decision_rules[p][1];
#         incoming_state = merge(Dict(Symbol("res_real[$(r)]") => r.currentvolume for r in res), Dict(Symbol("res_ind[$(r)]") => r.currentvolume for r in res)),
#         noise = (price = price_sample[rand(1:SCENARIO_COUNT)], inflow = 20),
#         controls_to_record = [:Qnom]
#     ).controls[:Qnom][r].out for r in res)
# end
# models = Dict(p => individual_solutions_anticipatory[p][2] for p in parts)
# Qnoms = ObtainNomination(decision_rules, parts, Dict(r => float(130 * 24 * 3) for r in res))
# print(Qnoms)