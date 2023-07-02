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
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation

filepath_Ljungan = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/Ljungan.json"
filepath_prices = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Data/Spot Prices/prices_df.csv"
filepath_results = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Results"
res, plants, parts = read_data(filepath_Ljungan)

function prepare_pricedata(filepath_prices)
    price_data = CSV.read(filepath_prices, DataFrame)
    price_data = coalesce.(price_data, 46.79)
    price_data.Sum = sum.(eachrow(price_data[:, 2:end]))/24
    rename!(price_data, :Column1 => :Date)
    return price_data
end

function create_Ω(price_data, inflow_scenarios, SCENARIO_COUNT::Int64, parts::Array{Participant}; quantile_bounds = 0.15, STAGE_COUNT = STAGE_COUNT)
    price_quantiles = quantile(price_data.Sum, range(quantile_bounds, 1 - quantile_bounds, length = SCENARIO_COUNT+1))
    price_subsets = Dict{Int64, DataFrame}()
    for i in 1:SCENARIO_COUNT
        price_subsets[i] = price_data[(price_data.Sum .>= price_quantiles[i]) .& (price_data.Sum .<= price_quantiles[i+1]), :]
    end
    price_scenarios = Dict{Int64, Vector{Any}}()
    for s in 1:SCENARIO_COUNT
        price_scenarios[s] = [collect(values(row)) for row in eachrow(select(price_subsets[s], Not([:Date, :Sum])))]
    end
    price_sample = Dict(s =>  [price_scenarios[scen][rand(1:length(price_scenarios[scen]))] for scen in 1:SCENARIO_COUNT] for s in 1:STAGE_COUNT)
    Ω = Dict(s => [(price = c, inflow = Q) for c in price_sample[s] for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
    P = Dict(s => [1/length(eachindex(Ω[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)
    Ω_scenario = Dict(scenario => Dict(s => [(price = Ω[s][scenario].price, inflow = Ω[s][scenario].inflow)] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)
    P_scenario = Dict(scenario => Dict(s => [1/length(eachindex(Ω_scenario[scenario][s])) for i in eachindex(Ω_scenario[scenario][s])] for s in 1:STAGE_COUNT) for scenario in 1:SCENARIO_COUNT)
    price_threshold::Dict{Participant, Dict{Reservoir, Float64}} = Dict(p => Dict(r => quantile(price_data.Sum,0.5) for r in res) for p in parts) 
    max_hourly_price = max([max([mean(s) for s in sample]...) for sample in values(price_sample)]...)
    return Ω, P, Ω_scenario, P_scenario, price_threshold, max_hourly_price, price_sample
end 

T = 24
STAGE_COUNT = 8
ITERATION_COUNT = 50
SCENARIO_COUNT = 10
BIG_M = 10e+4

price_data = prepare_pricedata(filepath_prices)
inflow_scenarios = [0]
Ω, P, Ω_scenario, P_scenario, mean_price, max_hourly_price, price_sample = create_Ω(price_data, inflow_scenarios, SCENARIO_COUNT, parts)
j = parts[2]
O,plants_O = OtherParticipant(parts, j, res)

INITIAL_RESERVOIR = Dict{Reservoir, Float64}(r => r.currentvolume for r in res)
INITIAL_INDIVIDUAL_RESERVOIR = Dict{Participant, Dict{Reservoir, Float64}}(p => Dict(r => p.individual_reservoir[r] for r in res) for p in parts)
Qref = Dict(r => mean(inflow_scenarios) for r in res)
nom = Dict(scenario => ShortTermOptimizationNoAnticipation(
        res,
        j,
        O,
        j.plants,
        ITERATION_COUNT,
        STAGE_COUNT,
        SCENARIO_COUNT,
        T,
        INITIAL_RESERVOIR,
        INITIAL_INDIVIDUAL_RESERVOIR[j],
        Ω_scenario[scenario],
        P_scenario[scenario],
        Qref,
        mean_price[j],
        price_sample;
        printlevel = 1
        )[3] for scenario in 1:SCENARIO_COUNT)

println(nom)
Ω_nom = Dict(s => [(price = price_sample[s][i], inflow = Q, nomination = nom[i][s]) for i in 1:length(price_sample[s]) for Q in inflow_scenarios] for s in 1:STAGE_COUNT)
P_nom = Dict(s => [1/length(eachindex(Ω_nom[s])) for i in eachindex(Ω[s])] for s in 1:STAGE_COUNT)

"""

"""
plot_safepath =" C://Users/lenna/OneDrive - NTNU/Master Thesis/Presentation VF/Images/othersNominations.png"
function plot_OtherNominations(nom)
    x = 1:length(nom[1])-1
    plot(x, [sum(nom[1][i][r] for r in res) for i in 1:length(nom[1])-1], legend=false, fmt=:png)
    for j in 2:length(nom)-1
        plot!(x, [sum(nom[j][i][r] for r in res) for i in 1:length(nom[1])-1], legend=false, fmt=:pmg)
    end
    xlabel!("Auction Period")
    ylabel!("Qnom")
    title!("Estimated Others' nomination")
    savefig(plot_safepath)
end


plot_OtherNominations(nom)

model, rules, nomination = ShortTermOptimizationAnticipation(
    res,
    j,
    O,
    plants_O,
    parts,
    j.plants,
    ITERATION_COUNT,
    STAGE_COUNT,
    SCENARIO_COUNT,
    T,
    INITIAL_RESERVOIR,
    INITIAL_INDIVIDUAL_RESERVOIR[j],
    Ω_nom,
    P_nom,
    nom,
    Qref,
    mean_price[j],
    price_sample;
    optimizer = CPLEX.Optimizer,
    printlevel = 1,
    BIG_M = 1e5,
    stall_bound = SDDP.BoundStalling(5, 1e-2))
nominations = []
rules = []
Qadjs = []
P_Overs = []
P_Swaps = []
for node in 1:STAGE_COUNT
    rule = SDDP.DecisionRule(model; node = node)
    solution = SDDP.evaluate(
        rule;
        incoming_state = merge(Dict(Symbol("res_real[$(r)]") => INITIAL_RESERVOIR[r] for r in res), Dict(Symbol("res_ind[$(r)]") =>INITIAL_INDIVIDUAL_RESERVOIR[j][r] for r in res)),
        noise = (price = price_sample[rand(1:STAGE_COUNT)][rand(1:SCENARIO_COUNT)], inflow = 0.0, nomination = nom[rand(1:SCENARIO_COUNT)]),
        controls_to_record = [:Qnom, :Qadj, :P_Over, :P_Swap]
    )
    Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
    Qadj = Dict(r => solution.controls[:Qadj] for r in res)
    P_Over = Dict(k => solution.controls[:P_Over] for k in plants_O)
    P_Swap = Dict(r => solution.controls[:P_Swap] for r in res)
    push!(rules, rule)
    push!(nominations, Qnom)
    push!(Qadjs, Qadj)
    push!(P_Overs, P_Over)
    push!(P_Swaps, P_Swap)
end

solution = SDDP.evaluate(
    rules[1];
    incoming_state = merge(Dict(Symbol("res_real[$(r)]") => INITIAL_RESERVOIR[r] for r in res), Dict(Symbol("res_ind[$(r)]") =>INITIAL_INDIVIDUAL_RESERVOIR[j][r] for r in res)),
    controls_to_record = [:Qnom]
)