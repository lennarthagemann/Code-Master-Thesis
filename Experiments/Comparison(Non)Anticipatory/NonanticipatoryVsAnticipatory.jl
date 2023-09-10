#=
-----------------------------------------------------------------------------------------------

We developed two different approaches to bidding in the water regulation 
We compare them side by side, while uncertainty reveals itself equaly.
How different are the profits? How does it affect each producer individually, whether other participants
use a certain strategy? Is any approach better than the other? Is it detrimental if everybody was 
nonanticipatory / anticipatory? 

-----------------------------------------------------------------------------------------------
=#


using JuMP
using CPLEX
using Distributions
using LinearAlgebra
using Statistics
using Dates
using DataFrames
using SDDP
using Plots
using CSV
using JSON
import Random: randperm
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end

includet(pwd() * "\\Water_Regulation\\WaterRegulation.jl")
using .WaterRegulation

const filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
const filepath_prices = pwd() * "\\Inflow Forecasting\\Data\\Spot Prices\\prices_df.csv"
const filepath_inflows = pwd() * "\\Inflow Forecasting\\Data\\Inflow\\Data from Flasjoen and Holmsjoen.csv"
const savepath_watervalue = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Inflow Forecasting\\WaterValue"

R, K, J = read_data(filepath_Ljungan)

const ColumnReservoir = Dict(r => r.dischargepoint * " Inflow" for r in R)
const scenario_count_inflows = 1
const scenario_count_prices = 10
const scenario_count_prices_medium = 3
const scenario_count_inflows_weekly = 3
const stage_count_short = 2
const stage_count_bidding = 2
const stage_count_medium = 52
const price_point_count = 5
const T = 24
const currentweek = 2
const iteration_count_short = 50
const iteration_count_Bidding = 100
const iteration_count_medium = 1000

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)

PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)
Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded = ReadMediumModel(savepath_watervalue, J, R, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
mu_up, mu_down = BalanceParameters(price_data)

Strategy_Combinations = Dict{Participant,String}[
    Dict(J[1] => "Nonanticipatory", J[2] => "Nonanticipatory", J[3] => "Nonanticipatory"),
    Dict(J[1] => "Nonanticipatory", J[2] => "Nonanticipatory", J[3] => "Anticipatory"),
    Dict(J[1] => "Nonanticipatory", J[2] => "Anticipatory", J[3] => "Nonanticipatory"),
    Dict(J[1] => "Nonanticipatory", J[2] => "Anticipatory", J[3] => "Anticipatory"),
    Dict(J[1] => "Anticipatory", J[2] => "Nonanticipatory", J[3] => "Nonanticipatory"),
    Dict(J[1] => "Anticipatory", J[2] => "Nonanticipatory", J[3] => "Anticipatory"),
    Dict(J[1] => "Anticipatory", J[2] => "Anticipatory", J[3] => "Nonanticipatory"),
    Dict(J[1] => "Anticipatory", J[2] => "Anticipatory", J[3] => "Anticipatory")
]

"""
function AnticipatoryVsNonanticipatory()

    Function in which we compare the Performance of all Strategies for each combination of Participants.
    We are interested in individual differences such as revenue and water usage from different reservoirs,
    but also system values, such as spillage and real water usage.
    Here we try to answer the question wether a single approach is more favorable than another.

"""
function AnticipatoryVsNonanticipatory(R::Vector{Reservoir},J::Vector{Participant}, mu_up::Float64, mu_down::Float64, inflow_data::DataFrame, price_data::DataFrame,
    MediumModel_j, MediumModel_O, currentweek::Int64, scenario_count_prices::Int64, iteration_count_bidding::Int64, Strategy::Dict{Participant, String}; ColumnReservoir = ColumnReservoir)
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
    Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
    cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
    WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
    WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)
    Qnoms_Bidding = Dict{Dict{Participant, String}, Dict{Participant, Dict{Reservoir, Float64}}}()
    Qnoms_Scheduling = Dict{Dict{Participant, String}, Dict{Participant, Dict{Reservoir, Float64}}}()
    Qadjs = Dict{Dict{Participant, String}, Dict{Reservoir, Float64}}()
    P_Swaps = Dict{Dict{Participant, String}, Dict{Participant, Dict{Reservoir, Float64}}}()
    z_ups = Dict{Dict{Participant, String}, Dict{Participant, Vector{Float64}} }()
    z_downs = Dict{Dict{Participant, String}, Dict{Participant, Vector{Float64}}}()
    Individual_Revenues = Dict{Dict{Participant, String}, Dict{Participant, Float64}}()
    l_reals = Dict{Dict{Participant, String}, Dict{Reservoir, Float64}}()
    l_inds = Dict{Dict{Participant, String}, Dict{Reservoir, Float64}}()
    for strat in Strategy_Combinations
        l_vf_ind = Dict(j => Dict(r => j.individualreservoir[r] for r in R) for j in J)
        HourlyBiddingCurves, Qnoms1, Ω1, PPoints = FirstLayerSimulation(R, strat, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek)
        price = Ω1[J[1]][stage_count_bidding][rand(1:scenario_count_prices)].price
        inflow = Dict(r => Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, scenario_count_inflows)[1][r][1] for r in R)
        Obligations = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
        Qadj1, _, P_Swap1, _, _, _ = water_regulation(Qnoms1, Qref, inflow, false)
        Qnoms2 = SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligations, price, price_data, inflow_data, Qref, cuts,  WaterCuts, iteration_count_short, mu_up, mu_down, T, stage_count_short)
        Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)
        z_up, z_down = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligations, mu_up, mu_down, T)
        Individual_Revenue = Final_Revenue(J, price, Obligations, z_up, z_down, mu_up, mu_down, T)
        for j in J
            for r in R
                l_vf_ind[j][r] = l_vf_ind[j][r] - Qnoms2[(participant = j, reservoir = r)] + Qref[r]
            end
        end
        
        Qnoms_Bidding[strat] = Qnoms1 
        Qnoms_Scheduling[strat] = Qnoms2
        Qadjs[strat] = Qadj2
        P_Swaps[strat] = P_Swap2
        z_ups[strat] = z_up
        z_downs[strat] = z_down
        Individual_Revenues[strat] = Individual_Revenue
        l_reals[strat] = Dict{Reservoir, Float64}(r => r.currentvolume for r in R)
        l_inds[strat] = l_vf_ind
    end

    return Qnoms, Qadjs, P_Swaps, z_ups, z_downs, Individual_Revenues, l_reals, l_inds
end

