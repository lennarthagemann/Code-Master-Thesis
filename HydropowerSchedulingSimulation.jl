# --------------------------------------------------- #
# In this file we simulate the different hydropower scheduling strategies over a time horizon.
# We obtain the nominations from optimization algorithms, feed them into the water regulation procedures to obtain the real results,
# and feed the real results into the optimization algorithms to obtain the next nominations.
# This is repeated over various time horizons and the results are compared.
# -------------------------1-Setup------------------- #
# Importing packages
# Read in reservoirs, participants, plants from json file.
# Define the current week, imply the season
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

"""
Reoptimize after Adjustment.
With Qadj and P_Swap as parameters, find an optimized schedule in the given timeframe T. 
"""
function get_future_season(season::String)::Array{String}
    if season == "Spring"
        return ["Summer"]
    elseif season == "Summer"
        return ["Autumn"]
    elseif season == "Autumn"
        return ["Winter"]
    else
        return ["Spring"]
    end
end

function get_season(month::Int64)
    if month in 3:4  
        return "Spring"
    elseif month in 5:8  
        return "Summer"
    elseif month in 9:11  
        return "Autumn"
    else  
        return "Winter"
    end
end

function get_weekend(weekday::Int64)
    if weekday in 1:5
        return "Weekday"
    else
        return "Weekend"
    end
end

function prepare_pricedata(filepath_prices)
    price_data = CSV.read(filepath_prices, DataFrame)
    price_data = coalesce.(price_data, 46.79)
    price_data.Sum = sum.(eachrow(price_data[:, 2:end]))/24
    rename!(price_data, :Column1 => :Date)
    price_data.season = get_season.(month.(price_data.Date))
    price_data.Weekday = dayofweek.(price_data.Date)
    price_data.Weekday = get_weekend.(price_data.Weekday)
    return price_data
end

function create_Ω(price_data, inflow_scenarios, weekday::String, season::String, SCENARIO_COUNT::Int64, parts::Array{Participant}; quantile_bounds = 0.1, STAGE_COUNT = STAGE_COUNT)
    price_quantiles = quantile(price_data.Sum, range(quantile_bounds, 1 - quantile_bounds, length = SCENARIO_COUNT+1))
    price_subsets = Dict{Int64, DataFrame}()
    for i in 1:SCENARIO_COUNT
        price_subsets[i] = price_data[(price_data.Sum .>= price_quantiles[i]) .& (price_data.Sum .<= price_quantiles[i+1]), :]
    end
    price_scenarios = Dict{Int64, Vector{Any}}()
    for s in 1:SCENARIO_COUNT
        price_scenarios[s] = [collect(values(row)) for row in eachrow(select(price_subsets[s], Not([:Date, :Weekday, :season, :Sum])))]
    end
    price_sample = Dict(s =>  [price_scenarios[scen][rand(1:length(price_scenarios[scen]))] for scen in 1:SCENARIO_COUNT] for s in 1:STAGE_COUNT)
    Ω = Dict(s => [(price = c, inflow = Q) for c in price_sample[s] for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
    P = Dict(s => [1/length(eachindex(Ω[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)
    Ω_scenario = Dict(scenario => Dict(s => [(price = Ω[s][scenario].price, inflow = Ω[s][scenario].inflow)] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)
    P_scenario = Dict(scenario => Dict(s => [1/length(eachindex(Ω_scenario[scenario][s])) for i in eachindex(Ω_scenario[scenario][s])] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)
    price_threshold::Dict{Participant, Dict{Reservoir, Float64}} = Dict(p => Dict(r => quantile(price_data.Sum,0.45) for r in res) for p in parts) 
    max_hourly_price = max([max([mean(s) for s in sample]...) for sample in values(price_sample)]...)
    return Ω, P, Ω_scenario, P_scenario, price_threshold, max_hourly_price, price_sample
end 

"""
Determine the price threshold of when to sell energy. This is determined by the percentage of water in the reservoir.
If it is only filled by 30% => producer power if price is higher that 70% of hourly prices.
"""

function ChooseStrategy(strat::String, p::Participant, r::Reservoir, Qnom, Qnom_ant)
    if strat == "Anticipatory"
        return Qnom_ant[(participant = p, reservoir = r)]
    else
        return Qnom[(participant = p, reservoir = r)]
    end
end

""" 
Hooray, time to optimize!
Obtain Nomination of first stage (as is the only deterministic one, and only realistic decision at the moment) for everyone by optimizing both strategies.
Return the models, rules and nominations.
"""
function OptimizationRound(
    INITIAL_RESERVOIR::Dict{Reservoir, Float64},
    INITIAL_INDIVIDUAL_RESERVOIR::Dict{Participant, Dict{Reservoir, Float64}},
    mean_hourly_price::Dict{Participant, Dict{Reservoir, Float64}};
    parts::Array{Participant} = parts,
    res::Array{Reservoir} = res,
    part_plants::Dict{Participant, Array{HydropowerPlant}} = part_plants,
    ITERATION_COUNT::Int64 = ITERATION_COUNT,
    STAGE_COUNT::Int64 = STAGE_COUNT,
    SCENARIO_COUNT::Int64 = SCENARIO_COUNT,
    T::Int64 = T,
    Ω = Ω,
    P = P,
    Qref = Qref,
    price_sample = price_sample,
    print_level = 0
    )
    nominations = Dict{Participant, Dict{Int64, Dict{Reservoir, Float64}}}(p => Dict{Int64, Dict{Reservoir, Float64}}() for p in parts)
    nominations_ant = Dict{Participant, Dict{Int64, Dict{Reservoir, Float64}}}(p => Dict{Int64, Dict{Reservoir, Float64}}() for p in parts)
    rules = Dict{Participant, Dict{Int64, Any}}(p => Dict() for p in parts)
    rules_ant = Dict{Participant, Dict{Int64, Any}}(p => Dict() for p in parts)
    models = Dict{Participant, Any}(p => Dict() for p in parts)
    models_ant = Dict{Participant, Any}(p => Dict() for p in parts)
    for p in parts
        local_model, local_rule, local_nom = ShortTermOptimizationNoAnticipation(
            res,
            p,
            OtherParticipant(parts, p, res)[1],
            part_plants[p],
            ITERATION_COUNT,
            STAGE_COUNT,
            SCENARIO_COUNT,
            T,
            INITIAL_RESERVOIR::Dict{Reservoir, Float64},
            INITIAL_INDIVIDUAL_RESERVOIR[p]::Dict{Reservoir, Float64},
            Ω,
            P,
            Qref,
            mean_hourly_price[p],
            price_sample;
            printlevel = print_level
        )
        models[p] = local_model
        rules[p] = Dict(i => local_rule[i] for i in eachindex(local_rule))
        nominations[p] = Dict(i => local_nom[i] for i in eachindex(local_nom))
    end

    for p in parts
        local_O, local_plants_O = OtherParticipant(parts, p ,res)

        local_nom = Dict(scenario => ShortTermOptimizationNoAnticipation(
        res,
        p,
        local_O,
        part_plants[p],
        ITERATION_COUNT,
        STAGE_COUNT,
        SCENARIO_COUNT,
        T,
        INITIAL_RESERVOIR::Dict{Reservoir, Float64},
        INITIAL_INDIVIDUAL_RESERVOIR[p]::Dict{Reservoir, Float64},
        Ω_scenario[scenario],
        P_scenario[scenario],
        Qref,
        mean_hourly_price[p],
        price_sample;
        printlevel = print_level
        )[3] for scenario in 1:SCENARIO_COUNT)

        Ω_nom_local = Dict(s => [(price = price_sample[s][i], inflow = Q, nomination = local_nom[i][s]) for i in 1:length(price_sample[s]) for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
        P_nom_local = Dict(s => [1/length(eachindex(Ω_nom_local[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)

        local_model_ant, local_rule_ant, local_nom_ant = ShortTermOptimizationAnticipation(
            res,
            p,
            local_O,
            local_plants_O,
            parts,
            part_plants[p],
            ITERATION_COUNT,
            STAGE_COUNT,
            SCENARIO_COUNT,
            T,
            INITIAL_RESERVOIR,
            INITIAL_INDIVIDUAL_RESERVOIR[p],
            Ω_nom_local,
            P_nom_local,
            local_nom,
            Qref,
            mean_hourly_price[p],
            price_sample;
            printlevel = print_level
        )
        models_ant[p] = local_model_ant
        rules_ant[p] = Dict(i => local_rule_ant[i] for i in eachindex(local_rule_ant))
        nominations_ant[p] = Dict(i => local_nom_ant[i] for i in eachindex(local_nom_ant))
    end

    Qnom = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    Qnom_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    for p in parts
        sol = SDDP.evaluate(
            rules[p][1];
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => INITIAL_RESERVOIR[r] for r in filter(r -> p.participationrate[r] > 0,res)), Dict(Symbol("res_ind[$(r)]") =>INITIAL_INDIVIDUAL_RESERVOIR[p][r] for r in  filter(r -> p.participationrate[r] > 0 ,res))),
            controls_to_record = [:Qnom]
        )
        sol_ant = SDDP.evaluate(
            rules_ant[p][1];
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => INITIAL_RESERVOIR[r] for r in filter(r -> p.participationrate[r] > 0,res)), Dict(Symbol("res_ind[$(r)]") =>INITIAL_INDIVIDUAL_RESERVOIR[p][r] for r in  filter(r -> p.participationrate[r] > 0 ,res))),
            controls_to_record = [:Qnom]
        )
        for r in res
            if p.participationrate[r] > 0
                Qnom[(participant = p, reservoir = r)] = sol.controls[:Qnom][r].out
                Qnom_ant[(participant = p, reservoir = r)] = sol_ant.controls[:Qnom][r].out
            else
                Qnom[(participant = p, reservoir = r)] = Qref[r]
                Qnom_ant[(participant = p, reservoir = r)] = Qref[r]
            end
        end
    end

    for node in 1:STAGE_COUNT
        for p in parts
            for r in res
                nominations[p][node][r] = round(nominations[p][node][r]; digits = 4)
                nominations_ant[p][node][r] = round(nominations_ant[p][node][r]; digits = 4)
            end
        end
    end
    return models, rules, Qnom, nominations, models_ant, rules_ant, Qnom_ant, nominations_ant
end


function SubsequentOptimizations(
    parts::Array{Participant},
    res::Array{Reservoir},
    rounds::Int64,
    INITIAL_RESERVOIR::Dict{Reservoir, Float64},
    INITIAL_INDIVIDUAL_RESERVOIR::Dict{Participant, Dict{Reservoir, Float64}},
    start_date::String,
    strategy,
    mean_hourly_price;
    Qref = Qref)

    Qnom_round::Array{Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}} = []
    Qadj_round::Array{Dict{Reservoir, Float64}} = []
    QadjTot_round::Array{Dict{Reservoir, Float64}} = []
    P_Swap_round::Array{Dict{Participant, Dict{Reservoir, Float64}}} = []
    POver_round::Array{Dict{Participant, Dict{HydropowerPlant, Float64}}} = []
    ΣPOver_round::Array{Dict{HydropowerPlant, Float64}} = [] 
    MaxEnergy_round::Array{Dict{HydropowerPlant, Float64}} = []
    date = Dates.Date(start_date, "yyyy-mm-dd")
    for round in 1:rounds
        current_date = date + Dates.Day(round - 1)
        current_weekday = dayofweek(current_date) <= 5 ? "Weekday" : "Weekend"
        current_season = get_season(Dates.month(current_date))
        println("Wir sind in Runde: ---------$(round)----------")
        Ω, P, Ω_scenario, P_scenario, mean_hourly_price, max_hourly_price, price_sample = create_Ω(price_data, inflow_scenarios, current_weekday, current_season, SCENARIO_COUNT, parts)
        if round == 1
            models, rules, Qnom, nominations_round, models_ant, rules_ant, Qnom_ant, nominations_ant = OptimizationRound(
            INITIAL_RESERVOIR,
            INITIAL_INDIVIDUAL_RESERVOIR,
            mean_hourly_price;
            Ω = Ω,
            P = P,
            price_sample = price_sample)
        else
            models, rules, Qnom, nominations_round, models_ant, rules_ant, Qnom_ant, nominations_ant = OptimizationRound(
            Dict(r => r.currentvolume for r in res),
            Dict(p => Dict(r => p.individual_reservoir[r] for r in res) for p in parts),
            mean_hourly_price;
            Ω = Ω,
            P = P,
            price_sample = price_sample)
        end
        Communicated_Nominations = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}(
        (participant = p, reservoir = r) => ChooseStrategy(strategy[(participant = p, reservoir = r)], p, r, Qnom, Qnom_ant) for p in parts for r in res)
        Qadj, QadjTot, P_Swap, POver, ΣPOver, MaxEnergy = water_regulation(Communicated_Nominations, Qref, T)
        append!(Qadj_round, [Qadj]) 
        append!(QadjTot_round,[QadjTot]) 
        append!(P_Swap_round, [P_Swap])
        append!(Qnom_round, [Communicated_Nominations])
        append!(POver_round, [POver])
        append!(ΣPOver_round, [ΣPOver])
        append!(MaxEnergy_round, [MaxEnergy])
    end
    return Qnom_round, Qadj_round, QadjTot_round, P_Swap_round, POver_round, ΣPOver_round, MaxEnergy_round
end

function save_results(filepath_results, filename, strategy, Qnom_round, Qadj_round, QadjTot_round, P_Swap_round, POver_round, ΣPOver_round, MaxEnergy_round)
    results = Dict("strategy" => strategy, "Qnom" => Qnom_round, "Qadj" => Qadj_round, "QadjTot" => QadjTot_round, "P_Swap" => P_Swap_round, "POver" => POver_round, "ΣPOver" => ΣPOver_round, "MaxEnergy" => MaxEnergy_round)
    json_results = JSON.json(results)
    open(filepath_results * "/" * filename, "w") do file  
        write(file, json_results)
    end
end

function read_results(filepath_results, filename, json_dictionary_parts, json_dictionary_res, json_dictionary_plants, parts, res, plants, rounds)
    json_data = JSON.parsefile(filepath_results * "/" * filename)
    Qnoms =  [Dict((participant = p, reservoir = r) => json_data["Qnom"][i]["(participant = $(json_dictionary_parts[p]), reservoir = $(json_dictionary_res[r]))"] for p in parts for r in res) for i in 1:rounds]
    Qadjs = [Dict(r => json_data["Qadj"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    QadjTots = [Dict(r => json_data["QadjTot"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    P_Swaps = [Dict(p => Dict(r => json_data["P_Swap"][i][json_dictionary_parts[p]][json_dictionary_res[r]] for r in res) for p in parts) for i in 1:rounds]
    POvers = [Dict(p => Dict(plant => json_data["POver"][i][json_dictionary_parts[parts[1]]][json_dictionary_plants[plants[3]]] for plant in plants) for p in parts) for i in 1:rounds]
    ΣPOvers = [Dict(plant => json_data["ΣPOver"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    MaxEnergys = [Dict(plant => json_data["MaxEnergy"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    return Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys
end

filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
filepath_prices = pwd() *  "\\Data\\Spot Prices\\prices_df.csv"
filepath_results = pwd() * "\\Results\\LambdaZero\\"
res, plants, parts = read_data(filepath_Ljungan)

T = 24
STAGE_COUNT = 7
ITERATION_COUNT = 10
SCENARIO_COUNT = 3

INITIAL_RESERVOIR = Dict{Reservoir, Float64}(r => r.currentvolume for r in res)
INITIAL_INDIVIDUAL_RESERVOIR = Dict{Participant, Dict{Reservoir, Float64}}(p => Dict(r => p.individual_reservoir[r] for r in res) for p in parts)

part_plants = Dict{Participant, Array{HydropowerPlant}}(p => p.plants for p in parts)

# ---------------------- 1 Probability Distributions ------------------- #
# 1. Get the big set of two dimensional uncertainty: prices and inflow
# 2. Define scenario-wise uncertainty sets: Deterministic, but different stage-wise
# 3.Calculate Uncertainty set for Anticipative model:
#   3.1 Get nominations for every possible price and inflow scenario
#   3.2 Make the uncertainty set 3-dimensional by adding the obtained nomination to the triple:
#       -> (price = ... , inflow = ... , nomination  = nom(price,inflow))

price_data = prepare_pricedata(filepath_prices)
inflow_scenarios = [0]
weekday = "Weekday"
season = "Summer"

Ω, P, Ω_scenario, P_scenario, mean_hourly_price, max_hourly_price = create_Ω(price_data, inflow_scenarios, weekday, season, SCENARIO_COUNT, parts)

Qref = Dict{Reservoir, Float64}(r => mean(inflow_scenarios) for r in res)


# -------------------------2-Optimization------------ #
# Use the chosen optimization technique for each participant
# Optimize based on the current parameters of the reservoir etc.
# Obtain the nominations for the next week
"""
Main Function for analysis. Sequentially obtain nominations, communicate nomination (with strategy) to VF.
Announce the adjusted flow and P_Swap, changed reservoir trajectories. (Subsequent Reoptimization to obtain real water flow)
Save the data from this round to return later.
Do a prescpecified amount of rounds. Return the data as Array with a length of rounds-
"""

strategy_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, String}(
    
    (participant = parts[1], reservoir = res[1]) => "Anticipatory",
    (participant = parts[1], reservoir = res[2]) => "Anticipatory",
    (participant = parts[2], reservoir = res[1]) => "Anticipatory",
    (participant = parts[2], reservoir = res[2]) => "Anticipatory",
    (participant = parts[3], reservoir = res[1]) => "Anticipatory",
    (participant = parts[3], reservoir = res[2]) => "Anticipatory")

strategy_no_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, String}(    
    (participant = parts[1], reservoir = res[1]) => "",
    (participant = parts[1], reservoir = res[2]) => "",
    (participant = parts[2], reservoir = res[1]) => "",
    (participant = parts[2], reservoir = res[2]) => "",
    (participant = parts[3], reservoir = res[1]) => "",
    (participant = parts[3], reservoir = res[2]) => "")

strategy_mixed1_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, String}(    
    (participant = parts[1], reservoir = res[1]) => "Anticipatory",
    (participant = parts[1], reservoir = res[2]) => "Anticipatory",
    (participant = parts[2], reservoir = res[1]) => "",
    (participant = parts[2], reservoir = res[2]) => "",
    (participant = parts[3], reservoir = res[1]) => "",
    (participant = parts[3], reservoir = res[2]) => "")

strategy_mixed2_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, String}(    
    (participant = parts[1], reservoir = res[1]) => "",
    (participant = parts[1], reservoir = res[2]) => "",
    (participant = parts[2], reservoir = res[1]) => "Anticipatory",
    (participant = parts[2], reservoir = res[2]) => "Anticipatory",
    (participant = parts[3], reservoir = res[1]) => "",
    (participant = parts[3], reservoir = res[2]) => "")

strategy_mixed3_ant = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, String}(    
    (participant = parts[1], reservoir = res[1]) => "",
    (participant = parts[1], reservoir = res[2]) => "",
    (participant = parts[2], reservoir = res[1]) => "",
    (participant = parts[2], reservoir = res[2]) => "",
    (participant = parts[3], reservoir = res[1]) => "Anticipatory",
    (participant = parts[3], reservoir = res[2]) => "Anticipatory")

    rounds = 7
start_date = "2020-04-03"

current_strat = strategy_mixed3_ant

Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys = SubsequentOptimizations(                             
        parts,
        res,                                       
        rounds,                                                                                                                                                                                                                                               
        INITIAL_RESERVOIR,                                                                                                            
        INITIAL_INDIVIDUAL_RESERVOIR,
        start_date,
        current_strat,                                                                                                 
        mean_hourly_price)
# ----  ---------------------3-Regulation-------------- #
# 1. Gather the nominations from every producer for every river
# 2. Organize the nominations, so that they can be fed into the subsequent calculations
# Feed the nominations into the water regulation procedure
# Calculate the adjusted flow, power swap and reduction -> Nominations as tuples (participant = ..., reservoir = ...) and Qref


# -------------------------4-Evaluate------------------ #
# How much Power was produced in total, how much belongs to each producer?
# How do the nominations differ from the expected other nomination?
# How much do the producers gain -> Fix a price scenario, and evaluate the revenue. (After Reoptimization)



# -------------------------5-Save-Results------------ #
# Save relevant results to track later and measure the performance of the different strategies.
# Is the anticipative approach smarter? When does it tend to have higher/lower profits?
# Is the system more competitive or less competitive?


"""
Save the results in .json format for further analysis.
"""
filename = "ResultsMixedStrat35.json"
save_results(filepath_results, filename, current_strat, Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys)


# -------------------------6-Plots--------------------- #
# Plot Nominations from the planning horizon, from either optimization algorithm.
# Dict(Participant ->[ Nominations ])
# Overlay Nominations with adjusted flow.