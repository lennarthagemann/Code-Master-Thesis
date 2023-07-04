using JSON
using Plots
using JuMP
using CPLEX
using CSV
using DataFrames
using Statistics
using Dates
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet(pwd() * "\\Water_Regulation\\WaterRegulation.jl")
using .WaterRegulation


filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
filepath_results = pwd() * "\\Results\\LambdaZero"
filenames_no_ant = ["ResultsNoAnticipation$(i).json" for i in 1:10]
filenames_ant = ["ResultsAnticipation$(i).json" for i in 1:10]
filenames_mixed_strat1 = ["ResultsMixedStrat1$(i).json" for i in 1:5]
filenames_mixed_strat2 = ["ResultsMixedStrat2$(i).json" for i in 1:5]
filenames_mixed_strat3 = ["ResultsMixedStrat3$(i).json" for i in 1:5]

res, plants, parts = read_data(filepath_Ljungan)
rounds = 7
T = 24
SCENARIO_COUNT=100
STAGE_COUNT = 1

json_dictionary_parts = Dict(
    parts[1] => "Name: Sydkraft\n",
    parts[2] => "Name: Fortum\n",
    parts[3] => "Name: Statkraft\n" )

json_dictionary_res = Dict(
    res[1] => "Reservoir with name: Flasjön\n",
    res[2] => "Reservoir with name: Holmsjön\n")

json_dictionary_plants = Dict(
    plants[1] => "A hydropowerplant with the name Flasjö and following information: \nName      : Flasjö\nReservoir : Flasjön\nEquivalent: 0.31\nSpill reference level: 0.58\n",
    plants[2] => "A hydropowerplant with the name Trangfors and following information: \nName      : Trangfors\nReservoir : Flasjön\nEquivalent: 0.72\nSpill reference level: 1.05\n",
    plants[3] => "A hydropowerplant with the name Rätan and following information: \nName      : Rätan\nReservoir : Flasjön\nEquivalent: 0.55\nSpill reference level: 1.15\n",
    plants[4] => "A hydropowerplant with the name Turinge and following information: \nName      : Turinge\nReservoir : Flasjön\nEquivalent: 0.19\nSpill reference level: 1.05\n",
    plants[5] => "A hydropowerplant with the name Bursnäs and following information: \nName      : Bursnäs\nReservoir : Flasjön\nEquivalent: 0.07\nSpill reference level: 1.05\n",
    plants[6] => "A hydropowerplant with the name Järnvägsforsen and following information: \nName      : Järnvägsforsen\nReservoir : Holmsjön\nEquivalent: 0.74\nSpill reference level: 1.45\n",
    plants[7] => "A hydropowerplant with the name Parteboda and following information: \nName      : Parteboda\nReservoir : Holmsjön\nEquivalent: 0.27\nSpill reference level: 1.4\n",
    plants[8] => "A hydropowerplant with the name Hermansboda and following information: \nName      : Hermansboda\nReservoir : Holmsjön\nEquivalent: 0.1\nSpill reference level: 1.4\n",
    plants[9] => "A hydropowerplant with the name Ljunga and following information: \nName      : Ljunga\nReservoir : Holmsjön\nEquivalent: 0.45\nSpill reference level: 1.45\n",
    plants[10] => "A hydropowerplant with the name Nederede and following information: \nName      : Nederede\nReservoir : Holmsjön\nEquivalent: 0.07\nSpill reference level: 2.4\n",
    plants[11] => "A hydropowerplant with the name Skallböle and following information: \nName      : Skallböle\nReservoir : Holmsjön\nEquivalent: 0.18\nSpill reference level: 2.65\n",
    plants[12] => "A hydropowerplant with the name Matfors and following information: \nName      : Matfors\nReservoir : Holmsjön\nEquivalent: 0.8\nSpill reference level: 2.5\n",
    plants[13] => "A hydropowerplant with the name Viforsen and following information: \nName      : Viforsen\nReservoir : Holmsjön\nEquivalent: 0.07\nSpill reference level: 1.6\n",)

function read_results(filepath_results, filename, json_dictionary_parts, json_dictionary_res, json_dictionary_plants, parts, res, plants, rounds)
    json_data = JSON.parsefile(filepath_results * "/" * filename)
    Qnoms =  [Dict((participant = p, reservoir = r) => json_data["Qnom"][i]["(participant = $(json_dictionary_parts[p]), reservoir = $(json_dictionary_res[r]))"] for p in parts for r in res) for i in 1:rounds]
    Qadjs = [Dict(r => json_data["Qadj"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    QadjTots = [Dict(r => json_data["QadjTot"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    P_Swaps = [Dict(p => Dict(r => json_data["P_Swap"][i][json_dictionary_parts[p]][json_dictionary_res[r]] for r in res) for p in parts) for i in 1:rounds]
    POvers = [Dict(p => Dict(plant => json_data["POver"][i][json_dictionary_parts[p]][json_dictionary_plants[plant]] for plant in plants) for p in parts) for i in 1:rounds]
    ΣPOvers = [Dict(plant => json_data["ΣPOver"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    MaxEnergys = [Dict(plant => json_data["MaxEnergy"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    return Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys
end

function get_results(filepath, filenames)
    Qnoms = []
    Qadjs = []
    QadjTots = []
    P_Swaps = []
    POvers = []
    ΣPOvers = []
    MaxEnergys = []
    for f in filenames
        Qnom, Qadj, QadjTot, P_Swap, POver, ΣPOver, MaxEnergy = read_results(filepath, f, json_dictionary_parts, json_dictionary_res, json_dictionary_plants, parts, res, plants, rounds)
        push!(Qnoms, Qnom)
        push!(Qadjs, Qadj)
        push!(QadjTots, QadjTot)
        push!(P_Swaps, P_Swap)
        push!(POvers, POver)
        push!(ΣPOvers, ΣPOver)
        push!(MaxEnergys, MaxEnergy)
    end
    return Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys
end
function create_Ω(price_data, inflow_scenarios, SCENARIO_COUNT::Int64, parts::Array{Participant}; quantile_bounds = 0.1, STAGE_COUNT = STAGE_COUNT)
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
    price_threshold::Dict{Participant, Dict{Reservoir, Float64}} = Dict(p => Dict(r => quantile(price_data.Sum,0.45) for r in res) for p in parts) 
    max_hourly_price = max([max([mean(s) for s in sample]...) for sample in values(price_sample)]...)
    return Ω, P, Ω_scenario, P_scenario, price_threshold, max_hourly_price, price_sample
end 
function prepare_pricedata(filepath_prices)
    price_data = CSV.read(filepath_prices, DataFrame)
    price_data = coalesce.(price_data, 46.79)
    price_data.Sum = sum.(eachrow(price_data[:, 2:end]))/24
    rename!(price_data, :Column1 => :Date)
    return price_data
end

filepath_prices = pwd() * "\\Data/Spot Prices\\prices_df.csv"
price_data = prepare_pricedata(filepath_prices)
Ω, P, Ω_scenario, P_scenario, mean_hourly_price, max_hourly_price, price_sample = create_Ω(price_data, [0.0], SCENARIO_COUNT, parts)

Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys = get_results(filepath_results, filenames_no_ant)
Qnoms_ant, Qadjs_ant, QadjTots_ant, P_Swaps_ant, POvers_ant, ΣPOvers_ant, MaxEnergys_ant = get_results(filepath_results, filenames_ant)
Qnoms_mix1, Qadjs_mix1, QadjTots_mix1, P_Swaps_mix1, POvers_mix1, ΣPOvers_mix1, MaxEnergys_mix1 = get_results(filepath_results, filenames_mixed_strat1)
Qnoms_mix2, Qadjs_mix2, QadjTots_mix2, P_Swaps_mix2, POvers_mix2, ΣPOvers_mix2, MaxEnergys_mix2 = get_results(filepath_results, filenames_mixed_strat2)
Qnoms_mix3, Qadjs_mix3, QadjTots_mix3, P_Swaps_mix3, POvers_mix3, ΣPOvers_mix3, MaxEnergys_mix3 = get_results(filepath_results, filenames_mixed_strat3)
ΣQnoms = Dict(p => [] for p in parts)
ΣQnoms_ant = Dict(p => [] for p in parts)
ΣQnoms_mix1 = Dict(p => [] for p in parts)
ΣQnoms_mix2 = Dict(p => [] for p in parts)
ΣQnoms_mix3 = Dict(p => [] for p in parts)
ΣQadjs = []
ΣQadjs_ant = []
ΣQadjs_mix1 = []
ΣQadjs_mix2 = []
ΣQadjs_mix3 = []
ΣP_Swaps = Dict(p => [] for p in parts)
ΣP_Swaps_ant = Dict(p => [] for p in parts)
ΣP_Swaps_mix1 = Dict(p => [] for p in parts)
ΣP_Swaps_mix2 = Dict(p => [] for p in parts)
ΣP_Swaps_mix3 = Dict(p => [] for p in parts)
#Overnomination and Ersmax should be done on average for each power plant-
SumΣPOvers = Dict(k => [] for k in plants)
SumΣPOvers_ant = Dict(k => [] for k in plants)
SumΣPOvers_mix1 = Dict(k => [] for k in plants)
SumΣPOvers_mix2 = Dict(k => [] for k in plants)
SumΣPOvers_mix3 = Dict(k => [] for k in plants)
ΣMaxEnergys = Dict(k => [] for k in plants)
ΣMaxEnergys_ant = Dict(k => [] for k in plants)
ΣMaxEnergys_mix1 = Dict(k => [] for k in plants)
ΣMaxEnergys_mix2 = Dict(k => [] for k in plants)
ΣMaxEnergys_mix3 = Dict(k => [] for k in plants)

for i in eachindex(Qnoms)
    for j in eachindex(Qnoms[i])
        for p in parts
            push!(ΣQnoms[p], sum(Qnoms[i][j][(participant = p, reservoir = r)] for r in res))
        end
    end
end
for i in eachindex(Qnoms_ant)
    for j in eachindex(Qnoms_ant[i])
        for p in parts
            push!(ΣQnoms_ant[p], sum(Qnoms_ant[i][j][(participant = p, reservoir = r)] for r in res))
        end
    end
end
for i in eachindex(Qnoms_mix1)
    for j in eachindex(Qnoms_mix1[i])
        for p in parts
            push!(ΣQnoms_mix1[p], sum(Qnoms_mix1[i][j][(participant = p, reservoir = r)] for r in res))
        end
    end
end
for i in eachindex(Qnoms_mix2)
    for j in eachindex(Qnoms_mix2[i])
        for p in parts
            push!(ΣQnoms_mix2[p], sum(Qnoms_mix2[i][j][(participant = p, reservoir = r)] for r in res))
        end
    end
end
for i in eachindex(Qnoms_mix3)
    for j in eachindex(Qnoms_mix3[i])
        for  p in parts
            push!(ΣQnoms_mix3[p], sum(Qnoms_mix3[i][j][(participant = p, reservoir = r)] for r in res))
        end
    end
end

for i in eachindex(Qadjs)
    for j in eachindex(Qadjs[i])
        push!(ΣQadjs, sum(Qadjs[i][j][r] for r in res))
    end
end
for i in eachindex(Qadjs_ant)
    for j in eachindex(Qadjs_ant[i])
        push!(ΣQadjs_ant, sum(Qadjs_ant[i][j][r] for r in res))
    end
end
for i in eachindex(Qadjs_mix1)
    for j in eachindex(Qadjs_mix1[i])
        push!(ΣQadjs_mix1, sum(Qadjs_mix1[i][j][r] for r in res))
    end
end
for i in eachindex(Qadjs_mix2)
    for j in eachindex(Qadjs_mix2[i])
        push!(ΣQadjs_mix2, sum(Qadjs_mix2[i][j][r] for r in res))
    end
end
for i in eachindex(Qadjs_mix3)
    for j in eachindex(Qadjs_mix3[i])
        push!(ΣQadjs_mix3, sum(Qadjs_mix3[i][j][r] for r in res))
    end
end

for i in eachindex(P_Swaps)
    for j in eachindex(P_Swaps[i])
        for p in parts
            push!(ΣP_Swaps[p], sum(P_Swaps[i][j][p][r] for r in res))
        end
    end
end
for i in eachindex(P_Swaps_ant)
    for j in eachindex(P_Swaps_ant[i])
        for p in parts
            push!(ΣP_Swaps_ant[p], sum(P_Swaps_ant[i][j][p][r] for r in res))
        end
    end
end
for i in eachindex(P_Swaps_mix1)
    for j in eachindex(P_Swaps_mix1[i])
        for p in parts
            push!(ΣP_Swaps_mix1[p], sum(P_Swaps_mix1[i][j][p][r] for r in res))
        end
    end
end
for i in eachindex(P_Swaps_mix2)
    for j in eachindex(P_Swaps_mix2[i])
        for p in parts
            push!(ΣP_Swaps_mix2[p], sum(P_Swaps_mix2[i][j][p][r] for r in res))
        end
    end
end
for i in eachindex(P_Swaps_mix3)
    for j in eachindex(P_Swaps_mix3[i])
        for p in parts
            push!(ΣP_Swaps_mix3[p], sum(P_Swaps_mix3[i][j][p][r] for r in res))
        end
    end
end

for i in eachindex(ΣPOvers)
    for j in eachindex(ΣPOvers[i])
        for p in plants
            push!(SumΣPOvers[p], ΣPOvers[i][j][p])
        end
    end
end
for i in eachindex(ΣPOvers_ant)
    for j in eachindex(ΣPOvers_ant[i])
        for p in plants
            push!(SumΣPOvers_ant[p], ΣPOvers_ant[i][j][p])
        end
    end
end
for i in eachindex(ΣPOvers_mix1)
    for j in eachindex(ΣPOvers_mix1[i])
        for p in plants
            push!(SumΣPOvers_mix1[p], ΣPOvers_mix1[i][j][p])
        end
    end
end
for i in eachindex(ΣPOvers_mix2)
    for j in eachindex(ΣPOvers_mix2[i])
        for p in plants
            push!(SumΣPOvers_mix2[p], ΣPOvers_mix2[i][j][p])
        end
    end
end
for i in eachindex(ΣPOvers_mix3)
    for j in eachindex(ΣPOvers_mix3[i])
        for p in plants
            push!(SumΣPOvers_mix3[p], ΣPOvers_mix3[i][j][p])
        end
    end
end

for i in eachindex(MaxEnergys)
    for j in eachindex(MaxEnergys[i])
        for p in plants
            push!(ΣMaxEnergys[p], sum(MaxEnergys[i][j][p]))
        end
    end
end
for i in eachindex(MaxEnergys_ant)
    for j in eachindex(MaxEnergys_ant[i])
        for p in plants
            push!(ΣMaxEnergys_ant[p], sum(MaxEnergys_ant[i][j][p]))
        end
    end
end
for i in eachindex(MaxEnergys_mix1)
    for j in eachindex(MaxEnergys_mix1[i])
        for p in plants
            push!(ΣMaxEnergys_mix1[p], sum(MaxEnergys_mix1[i][j][p]))
        end
    end
end
for i in eachindex(MaxEnergys_mix2)
    for j in eachindex(MaxEnergys_mix2[i])
        for p in plants
            push!(ΣMaxEnergys_mix2[p], sum(MaxEnergys_mix2[i][j][p]))
        end
    end
end
for i in eachindex(MaxEnergys_mix3)
    for j in eachindex(MaxEnergys_mix3[i])
        for p in plants
            push!(ΣMaxEnergys_mix3[p], sum(MaxEnergys_mix3[i][j][p]))
        end
    end
end

println("mean amount of nominations (Nonanticipatory): ", Dict(k => [mean(v)] for (k,v) in ΣQnoms))
println("mean amount of nominations (Anticipatory): ", Dict(k => [mean(v)] for (k,v) in ΣQnoms_ant))
println("mean amount of nominations (Mixed) Participant 1: ", Dict(k => [mean(v)] for (k,v) in ΣQnoms_mix1))
println("mean amount of nominations (Mixed) Participant 2: ", Dict(k => [mean(v)] for (k,v) in ΣQnoms_mix2))
println("mean amount of nominations (Mixed) Participant 3: ", Dict(k => [mean(v)] for (k,v) in ΣQnoms_mix3))
println("________________________________")
println("mean amount of adjusted Flow (Nonanticipatory): ", mean(ΣQadjs))
println("mean amount of adjusted Flow (Anticipatory): ", mean(ΣQadjs_ant))
println("mean amount of adjusted Flow (Mixed) Participant 1: ", mean(ΣQadjs_mix1))
println("mean amount of adjusted Flow (Mixed) Participant 2: ", mean(ΣQadjs_mix2))
println("mean amount of adjusted Flow (Mixed) Participant 3: ", mean(ΣQadjs_mix3))
println("________________________________")
println("mean Power Swap(Nonanticipatory): ", Dict(k => [mean(v)] for (k,v) in ΣP_Swaps ))
println("mean Power Swap (Anticipatory): ", Dict(k => [mean(v)] for (k,v) in ΣP_Swaps_ant ))
println("mean Power Swap (Mixed) Participant 1: ", Dict(k => [mean(v)] for (k,v) in ΣP_Swaps_mix1 ))
println("mean Power Swap (Mixed) Participant 2: ", Dict(k => [mean(v)] for (k,v) in ΣP_Swaps_mix2 ))
println("mean Power Swap (Mixed) Participant 3: ", Dict(k => [mean(v)] for (k,v) in ΣP_Swaps_mix3 ))
println("________________________________")
println("mean amount of total overnomination: ", Dict(k => mean(v) for (k,v) in SumΣPOvers))
println("mean amount of total overnomination (Anticipatory): ", Dict(k => mean(v) for (k,v) in SumΣPOvers_ant))
println("mean amount of total overnomination (Mixed) Participant 1: ", Dict(k => mean(v) for (k,v) in SumΣPOvers_mix1))
println("mean amount of total overnomination (Mixed) Participant 2: ", Dict(k => mean(v) for (k,v) in SumΣPOvers_mix2))
println("mean amount of total overnomination (Mixed) Participant 3: ", Dict(k => mean(v) for (k,v) in SumΣPOvers_mix3))
println("________________________________")
println("mean amount of lost energy: ", Dict(k => mean(v) for (k,v) in ΣMaxEnergys))
println("mean amount of lost energy (Anticipatory): ", Dict(k => mean(v) for (k,v) in ΣMaxEnergys_ant))
println("mean amount of lost energy (Mixed) Participant 1: ", Dict(k => mean(v) for (k,v) in ΣMaxEnergys_mix1))
println("mean amount of lost energy (Mixed) Participant 2: ", Dict(k => mean(v) for (k,v) in ΣMaxEnergys_mix2))
println("mean amount of lost energy (Mixed) Participant 3: ", Dict(k => mean(v) for (k,v) in ΣMaxEnergys_mix3))

function share_of_perfect_value(c::Vector{Float64}, parts, Qnoms, Qadjs, P_Swaps)
    share = Dict(p => Dict(j => Dict( i => 0.0 for i in 1:rounds) for j in 1:length(Qadjs)) for p in parts)
    for j in eachindex(Qadjs)
        for p in parts
            for i in 1:rounds
                Qeff_test_real, real_value = OptimizationAfterAdjustment(res, p.plants, Qadjs[j][i], P_Swaps[j][i][p], c, T)
                Qeff_test_theoretical, theoretical_value = OptimizationAfterAdjustment(res, p.plants, Dict(r => Qnoms[j][i][(participant = p, reservoir = r)] for r in res), Dict(r => 0.0 for r in res), c, T)
                if isapprox(real_value,0.0,atol=1e-6)
                    share[p][j][i] = 1.0
                else
                    share[p][j][i] = real_value/theoretical_value
                end
            end
        end
    end
    means = Dict(p => 0.0 for p in parts)
    for (p, inner_dict) in share
        sum_values = 0.0
        count_values = 0
        for (_, shares) in inner_dict
            sum_values += sum(values(shares))
            count_values += length(shares)
        end
        means[p] = sum_values / count_values
    end
    return means
end

function total_share(prices::Vector{Vector{Float64}}, parts,  Qnoms, Qadjs, P_Swaps)
    tot_share = Dict(p => 0.0 for  p in parts)
    for price in prices
        share = share_of_perfect_value(price, parts, Qnoms, Qadjs, P_Swaps)
        for p in parts
            tot_share[p] += share[p]
        end
    end
    for p in parts
        tot_share[p] = tot_share[p]/length(prices)
    end
    return tot_share
end
c= [float(i) for i in 1:T]

# means_no_ant = total_share(price_sample[1],parts, Qnoms, Qadjs, P_Swaps)
# means_ant = total_share(price_sample[1], parts, Qnoms_ant, Qadjs_ant, P_Swaps_ant)
# means_mix1 = total_share(price_sample[1], parts, Qnoms_mix1, Qadjs_mix1, P_Swaps_mix1)
# means_mix2 = total_share(price_sample[1], parts, Qnoms_mix2, Qadjs_mix2, P_Swaps_mix2)
# means_mix3 = total_share(price_sample[1], parts, Qnoms_mix3, Qadjs_mix3, P_Swaps_mix3)
