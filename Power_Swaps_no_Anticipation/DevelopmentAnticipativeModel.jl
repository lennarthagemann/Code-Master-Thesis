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
# Hydopower Data
filepath_systemA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/Ljungan.json"
res, plants, parts = read_data(filepath_systemA)
j = parts[3]
part_plants = Dict{Participant, Array{HydropowerPlant}}(p => p.plants for p in parts) 
plants_j = part_plants[j]
O, plants_O = OtherParticipant(parts, j, res)


# Parameters 
stall_bound = SDDP.BoundStalling(5, 1e-4)
T = 24
STAGE_COUNT = 8
ITERATION_COUNT = 50
SCENARIO_COUNT = 3
BIG_M = 10e+4
INITIAL_RESERVOIR = Dict(r => r.currentvolume for r in res)
INITIAL_INDIVIDUAL_RESERVOIR = Dict(r => r.currentvolume for r in res)

filepath_prices = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Data/Spot Prices/prices_df.csv"
price_data = CSV.read(filepath_prices, DataFrame)
price_data = coalesce.(price_data, 46.79)
# Get Array of every row in price_data (excluding Data Column1)
price_scenarios = []
for row in eachrow(select(price_data, Not(:Column1)))
    push!(price_scenarios, collect(values(row)))
end
price_sample = Dict(s => price_scenarios[rand(1:length(price_scenarios), SCENARIO_COUNT)] for s in 1:STAGE_COUNT)
max_hourly_price = max([max([mean(s) for s in sample]...) for sample in values(price_sample)]...)
mean_hourly_price = mean(mean(mean(mean(sample, dims=2)) for sample in values(price_sample)))
inflow_scenarios = [0.1]
Qref = Dict{Reservoir, Float64}(r => mean(inflow_scenarios) for r in res)
Ω = Dict(s => [(price = c, inflow = Q) for c in price_sample[s] for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
P = Dict(s => [1/length(eachindex(Ω[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)
#2
Ω_scenario = Dict(scenario => Dict(s => [(price = Ω[s][scenario].price, inflow = Ω[s][scenario].inflow)] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)
P_scenario = Dict(scenario => Dict(s => [1/length(eachindex(Ω_scenario[scenario][s])) for i in eachindex(Ω_scenario[scenario][s])] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)

#3.1
nom = Dict(scenario => ShortTermOptimizationNoAnticipation(
    res,
    j,
    OtherParticipant(parts, j, res)[1],
    part_plants[j],
    ITERATION_COUNT,
    STAGE_COUNT,
    SCENARIO_COUNT,
    T,
    INITIAL_RESERVOIR::Dict{Reservoir, Float64},
    Dict(r => r.currentvolume for r in res),
    Ω_scenario[scenario],
    P_scenario[scenario],
    Qref,
    mean_hourly_price,
    price_sample;
    printlevel = false
    )[3] for scenario in 1:SCENARIO_COUNT)
#3.2
# (Link price scenario from each stage with the nomination from a scenario at a certain stage (Mind Blown))
Ω_nom = Dict(s => [(price = price_sample[s][i], inflow = Q, nomination = nom[i][s]) for i in 1:length(price_sample[s]) for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
P_nom = Dict(s => [1/length(eachindex(Ω_nom[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)




"""
This function is a deterministic one-stage model of the water regulation problem.
Its purpose is to provide a guess for the other participants' decisions in the current stage.
It is parametrized by the scenario prize and inflow, as well as the current reservoir filling.
"""
function DeterministicNoAnticipation(
    res::Array{Reservoir},
    plants::Array{HydropowerPlant},
    price::Vector{Float64},
    Qinflow::Dict{Reservoir, Float64},
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Qref::Dict{Reservoir, Float64}
    )
    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_initial[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_initial[r])
        # Control Variables
        @variable(subproblem, Qnom[r = res])
        @variable(subproblem, Qeff[k = plants, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        
        for r in res
            @constraint(subproblem, Qnom[r] >= - sum(Qnom[d] for d in filter(x -> x != r, find_us_reservoir(r))))
            @constraint(subproblem, sum(Qreal[r,t] for t in 1:T) == T * Qnom[r])
            @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r] - Qinflow[r]))
            @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r] - Qref[r]))
            for k in filter(k -> k.reservoir == r, plants)
                @constraint(subproblem, sum(Qnom[d] for d in find_us_reservoir(r)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
            end
            @constraint(subproblem, Qnom[r] <= Qref[r] + BIG_M * (1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
        end
        for k in plants
            for t in 1:T
                @constraint(subproblem, Qeff[k,t] <= sum(Qreal[d, t] for d in find_us_reservoir(k.reservoir)))
                @constraint(subproblem, Qeff[k,t] <= k.spill_reference_level)
            end
        end

        @stageobjective(subproblem, -(sum(price[t] * Qeff[k,t] * k.equivalent for t in 1:T for k in plants)
      + sum((res_real[r].out - res_real[r].in) * sum(k.equivalent for k in filter(x -> x.reservoir in find_ds_reservoirs(r), plants)) for r in res) * mean_hourly_price))
        return subproblem
    end

    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = 1,
        sense = :Min,
        lower_bound = -sum(res_real_initial[r] * max_hourly_price for r in res),
        optimizer = CPLEX.Optimizer
    )
    SDDP.train(model; iteration_limit = 5, print_level = 0)
    rule = SDDP.DecisionRule(model; node = 1)
    solution = SDDP.evaluate(
        rule;
        incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
        controls_to_record = [:Qnom, :Qreal]
    )
    Qnom_O = Dict(r => solution.controls[:Qnom][r] for r in res)
    return Qnom_O
end

function ShortTermOptimizationAnticipationTest(
    res::Array{Reservoir},
    parts::Array{Participant},
    j::Participant,
    plants_j::Array{HydropowerPlant},
    ITERATION_COUNT::Int64,
    stage_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Qref = Qref;
    Ω = Ω_nom,
    P = P_nom,
    optimizer = CPLEX.Optimizer
    )
    O, plants_O = OtherParticipant(parts, j, res)
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
        @variable(subproblem, Qnom_O[r = res])
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
                    try
                        JuMP.fix(Qnom_O[r], ω.nomination[r], force=true)
                    catch
                        println("Die Datentypen stimmen nicht überein, die Nomination muss eine einzelne Zahl sein.")
                        for r in res
                            println(ω.nomination[r])
                        end
                    end
                end
            end 
            for r in res
                # Transition Function
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qadj[r]- Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                # Constraints
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
                @constraint(subproblem, Qadj[r] == (Qnom_O[r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
                @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
                @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
                for k in plants_O
                    @constraint(subproblem, Qnom_change[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
                end
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            for k in plants_O
                @constraint(subproblem, P_Over[k] >=  (Qadj[k.reservoir] - k.spill_reference_level) * k.equivalent)
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
            + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price
            + sum((1 - (j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_ind[r].out - res_ind[r].in)  * j.participationrate[r] for r in res) * mean_hourly_price))
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
    SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [SDDP.BoundStalling(2, 1e-4)], print_level = 1)
    # obtain decision rule in all steps, as well as nominations for the reservoir situation.
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 20, nomination = Dict(r => nom[rand(1:SCENARIO_COUNT)] for r in res)),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        push!(rules, rule)
        push!(nominations, Qnom)
    end
    return model, rules, nominations
end
"""
Diese Funktion ist nur logisch für ein System mit einem einzelnen Reservoir!
(Ich hasse diese blöde Kaskade, es könnte alles so schön sein)
"""

function ShortTermOptimizationNoAnticipationTest(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    ITERATION_COUNT::Int64,
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
                for k in plants
                    @constraint(subproblem, Qnom_change[r] <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[r]))
                end
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            for t in 1:T
                for k in plants_j
                    @constraint(subproblem, Qeff[k, t] <= Qreal[k.reservoir, t])
                    @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
                end
            end
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
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
    SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound])
    # obtain decision rule in first step
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 20),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        append!(rules, [rule])
        append!(nominations, [Qnom])
    end
    return model, rules, nominations
end

"""
For cascaded reservoirs, certain restrictions have to be adapted.
-> For spillage / effective flow the total flow has to be considered (in contrast to the local nominations)
"""
function ShortTermOptimizationNoAnticipationCascadeTest(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    ITERATION_COUNT::Int64,
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
        SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound])
    else
        SDDP.train(model; iteration_limit = ITERATION_COUNT, print_level = 1, stopping_rules = [stall_bound])
    end
    # obtain decision rule in first step
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 20),
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

function ShortTermOptimizationAnticipationCascadeTest(
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
    mean_price,
    Qref = Qref;
    optimizer = CPLEX.Optimizer
    )
    println(O.participationrate)
    println(j.participationrate)
    function subproblem_builder(subproblem::Model, node::Int)
        # State Variables
        @variable(subproblem, 0 <= res_real[r = res] <= r.maxvolume, SDDP.State, initial_value = res_real_initial[r])
        @variable(subproblem, res_ind[r = res], SDDP.State, initial_value = res_ind_initial[r])
        @variable(subproblem, Qnom[r = res], SDDP.State, initial_value = 0)
        # Control Variables
        @variable(subproblem, Qnom_change[r = res] >= 0)
        @variable(subproblem, Qeff[k = plants_j, t = 1:T] >= 0)
        @variable(subproblem, Qreal[r = res, t = 1:T] >= 0)
        @variable(subproblem, Qadj[r = res] >= 0)
        @variable(subproblem, Qmax[k = plants_O] >= 0)
        @variable(subproblem, P_Swap[r = res])
        @variable(subproblem, P_Over[k = plants_O] >= 0)
        @variable(subproblem, BALANCE_INDICATOR[r = res], Bin)
        # Random Variables
        @variable(subproblem, c[t = 1:T])
        @variable(subproblem, Qinflow[r = res] >= 0)
        @variable(subproblem, Qnom_O[r = res])
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
                    try
                        JuMP.fix(Qnom_O[r], ω.nomination[r], force=true)
                    catch
                        println("Die Datentypen stimmen nicht überein, die Nomination muss eine einzelne Zahl sein.")
                        for r in res
                            println(ω.nomination[r])
                        end
                    end
                end
            end 
            for r in res
                # Transition Function
                @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qadj[r]- Qinflow[r]))
                @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
                @constraint(subproblem, Qnom[r].out == Qnom_change[r])
                # Constraints
                @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
                @constraint(subproblem, Qadj[r] == (Qnom_O[r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
                @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
                @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
                @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
                @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            end
            for k in plants_O
                @constraint(subproblem, P_Over[k] >=  (Qadj[k.reservoir] - k.spill_reference_level) * k.equivalent)
                @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            end
            for t in 1:T
                for k in plants_j
                    @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                    @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
                end
            end
            # Objective Function
            @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
            + sum(c[t] * P_Swap[r] for t in 1:T for r in res)
            + sum((1 - (j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] for r in res) * mean_price
            + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r]))* (res_ind[r].out - res_ind[r].in)  * j.participationrate[r] for r in res) * mean_price)/1e3)
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
    SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound], risk_measure = SDDP.WorstCase() ,print_level = 1)
    # obtain decision rule in all steps, as well as nominations for the reservoir situation.
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 0.1, nomination = Dict(r => nom[rand(1:SCENARIO_COUNT)] for r in res)),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        push!(rules, rule)
        push!(nominations, Qnom)
    end
    return model, rules, nominations
end


OptimizedModelsParticipantsNoAnticipationCascade = Dict(p =>ShortTermOptimizationNoAnticipation(
    res,
    j,
    OtherParticipant(parts, j, res)[1],
    part_plants[j],
    ITERATION_COUNT,
    STAGE_COUNT,
    SCENARIO_COUNT,
    T,
    INITIAL_RESERVOIR::Dict{Reservoir, Float64},
    Dict(r => r.currentvolume for r in res),
    Ω,
    P,
    Qref,
    mean_hourly_price,
    price_sample;
    printlevel = false
)
for p in parts)

# OptimizedModelsParticipantsNoAnticipation = Dict(p =>ShortTermOptimizationNoAnticipationTest(
#     res,
#     p,
#     OtherParticipant(parts, p, res)[1],
#     filter(x -> x in p.plants, plants),
#     ITERATION_COUNT,
#     stage_count,
#     T,
#     INITIAL_RESERVOIR,
#     INITIAL_INDIVIDUAL_RESERVOIR,
#     true)
# for p in parts)

OptimizedModelsParticipantsAnticipationCascade = Dict(p => ShortTermOptimizationAnticipationCascadeTest(
    res,
    p,
    O,
    plants_O,
    parts,
    filter(x -> x in p.plants, plants),
    ITERATION_COUNT,
    STAGE_COUNT,
    SCENARIO_COUNT,
    T,
    INITIAL_RESERVOIR,
    INITIAL_INDIVIDUAL_RESERVOIR,
    Ω_nom,
    P_nom,
    mean_hourly_price)
for p in parts)

# OptimizedModelsParticipantsAnticipation = Dict(p =>ShortTermOptimizationAnticipationTest(
#     res,
#     parts,
#     p,
#     filter(x -> x in p.plants, plants),
#     ITERATION_COUNT,
#     STAGE_COUNT,
#     T,
#     INITIAL_RESERVOIR,
#     INITIAL_INDIVIDUAL_RESERVOIR)
# for p in parts)

# Get the decision rules for every stage, save them in rules Array

sims_no_anticipation_cascade = SDDP.simulate(
    OptimizedModelsParticipantsNoAnticipationCascade[parts[1]][1],
    100,
    [:Qnom, :Qreal, :Qeff, :res_real, :res_ind]
)

# sims_no_anticipation = SDDP.simulate(
#     OptimizedModelsParticipantsNoAnticipation[parts[1]][1],
#     100,
#     [:Qnom, :Qreal, :Qeff, :res_real, :res_ind]
# )

sims_anticipation_cascade = SDDP.simulate(
    OptimizedModelsParticipantsAnticipationCascade[parts[1]][1],
    100,
    [:Qnom_change, :Qreal, :Qeff, :res_real, :res_ind, :Qadj, :Qnom, :P_Swap, :P_Over, :Qnom_O]
)

# sims_anticipation = SDDP.simulate(
#     OptimizedModelsParticipantsAnticipation[parts[1]][1],
#     100,
#     [:Qnom_change, :Qreal, :Qeff, :res_real, :res_ind, :Qadj, :Qnom, :Qnom_O, :P_Swap, :P_Over]
# )



plt_no_anticipation_cascade = SDDP.SpaghettiPlot(sims_no_anticipation_cascade)
for r in res
    SDDP.add_spaghetti(plt_no_anticipation_cascade; title= "Reservoir volume - $(r.dischargepoint)") do data
        return data[:res_real][r].out
    end
    SDDP.add_spaghetti(plt_no_anticipation_cascade; title= "Individual Reservoir volume - $(r.dischargepoint) - $(j.name)") do data
        return data[:res_ind][r].out
    end
    SDDP.add_spaghetti(plt_no_anticipation_cascade; title = "Nomination - $(r.dischargepoint) - $(j.name)") do data
        return data[:Qnom][r].in
    end
end

SDDP.plot(plt_no_anticipation_cascade, "spaghetti_plot_not_anticipative.html")

# plt_no_anticipation = SDDP.SpaghettiPlot(sims_no_anticipation)
# for r in res
#     SDDP.add_spaghetti(plt_no_anticipation; title= "Reservoir volume - $(r.dischargepoint)") do data
#         return data[:res_real][r].out
#     end
#     SDDP.add_spaghetti(plt_no_anticipation; title= "Individual Reservoir volume - $(r.dischargepoint) - $(j.name)") do data
#         return data[:res_ind][r].out
#     end
#     SDDP.add_spaghetti(plt_no_anticipation; title = "Nomination - $(r.dischargepoint) - $(j.name)") do data
#         return data[:Qnom][r].in
#     end
# end

# SDDP.plot(plt_no_anticipation, "spaghetti_plot_not_anticipative.html")

plt_anticipation_cascade = SDDP.SpaghettiPlot(sims_anticipation_cascade)
for r in res
    SDDP.add_spaghetti(plt_anticipation_cascade; title= "Reservoir volume - $(r.dischargepoint)") do data
        return data[:res_real][r].out
    end
    SDDP.add_spaghetti(plt_anticipation_cascade; title= "Individual Reservoir volume - $(r.dischargepoint) - $(j.name)") do data
        return data[:res_ind][r].out
    end
    SDDP.add_spaghetti(plt_anticipation_cascade; title = "Nomination - $(r.dischargepoint) - $(j.name)") do data
        return data[:Qnom][r].in
    end
    SDDP.add_spaghetti(plt_anticipation_cascade; title = "Nomination - $(r.dischargepoint) - $(O.name)") do data
        return data[:Qnom_O][r]
    end
    SDDP.add_spaghetti(plt_anticipation_cascade; title= "Adjusted Flow - $(r.dischargepoint)") do data
        return data[:Qadj][r]
    end
    SDDP.add_spaghetti(plt_anticipation_cascade; title= "Power Swap- $(r.dischargepoint)") do data
        return data[:P_Swap][r]
    end
end
for k in plants_O
    SDDP.add_spaghetti(plt_anticipation_cascade; title= "Overnomination - $(k.reservoir.dischargepoint)") do data
        return data[:P_Over][k]
    end
end
SDDP.plot(plt_anticipation_cascade, "spaghetti_plot_anticipative.html")

# plt_anticipation = SDDP.SpaghettiPlot(sims_anticipation)
# for r in res
#     SDDP.add_spaghetti(plt_anticipation; title= "Reservoir volume - $(r.dischargepoint)") do data
#         return data[:res_real][r].out
#     end
#     SDDP.add_spaghetti(plt_anticipation; title= "Individual Reservoir volume - $(r.dischargepoint) - $(j.name)") do data
#         return data[:res_ind][r].out
#     end
#     SDDP.add_spaghetti(plt_anticipation; title = "Nomination - $(r.dischargepoint) - $(j.name)") do data
#         return data[:Qnom][r].in
#     end
#     SDDP.add_spaghetti(plt_anticipation; title = "Nomination - $(r.dischargepoint) - $(j.name)") do data
#         return data[:Qnom_O][r]
#     end
#     SDDP.add_spaghetti(plt_anticipation; title= "Adjusted Flow - $(r.dischargepoint)") do data
#         return data[:Qadj][r]
#     end
#     SDDP.add_spaghetti(plt_anticipation; title= "Power Swap- $(r.dischargepoint)") do data
#         return data[:P_Swap][r]
#     end
# end
# SDDP.plot(plt_anticipation, "spaghetti_plot_anticipative.html")

# V = SDDP.ValueFunction(model_anticipative, node = 1)
# V_anticipative = SDDP.evaluate(
    #     V, merge(Dict(Symbol("res_real[$(r)]") => INITIAL_RESERVOIR + 1000.0 for r in res),
    #              Dict(Symbol("res_ind[$(r)]") => INITIAL_INDIVIDUAL_RESERVOIR for r in res)))
    
    # Values_anticipative = [SDDP.evaluate(
#     V, merge(Dict(Symbol("res_real[$(r)]") => i  for r in res),
#              Dict(Symbol("res_ind[$(r)]") => i for r in res)))[1] for i in 1:Int(res[1].maxvolume)]

# plt = plot(1:Int(res[1].maxvolume), Values_anticipative, title = "Trial values of the value function")
# display(plt)