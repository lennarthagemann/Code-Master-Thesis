
#=
---------------------------------------------------------------------------------------------

If we consider the actions that happen over the course of 48 hours, the sequential decisions
can be modelled by the connection of the optimization models developed in this repository.

    - How much revenue is generated for each Participant after one simulation?
    - How much revenue is generated for each Participant after successive simulations? 
        For example one week/ one month
    - How do the results differ if we use the anticipatory model / nonanticipatory model?

---------------------------------------------------------------------------------------------
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
using UUIDs
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
NameToParticipant = Dict{String, Participant}(j.name => j for j in J)

PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)
Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded = ReadMediumModel(savepath_watervalue, J, R, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
MediumModelSingle = ReadMediumModelSingle(savepath_watervalue, R, K, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
Strategy = Dict(j => "Anticipatory" for j in J)
mu_up, mu_down = BalanceParameters(price_data)


"""
    Test_Price_Points()

    Create Price Points fitting to the scenarios generated.
    It is important to have a bidding curve fitting to all scenarios generated. Between every price point at least one scenario has to situated.
    Otherwise, the agent is indifferent to the volume to nominate at that point and a realized price would be detrimental.
    For every hour

"""
function Test_Price_Points(Ω, scenario_count_prices, T, max_point)::Vector{Vector{Float64}}
    if scenario_count_prices == 1
        PPoints = [[0.0, max_point] for t in 1:T]
    else
        PPoints = []
        for t in 1:T
            hourly_prices = sort(scen.price[t] for scen in Ω[2])
            temp = [(hourly_prices[i] + hourly_prices[i+1])/2 for i in 1:length(hourly_prices)-1]
            push!(PPoints, [0.0, temp..., max_point])
        end
    end
    return PPoints
end


"""
    SingleOwnerSimulation(R, K, price_data, inflow_data, cuts, WaterCuts, PPoints, iteration_count, mu_up, mu_down, T, stage_count_bidding)

    Simulate the Bidding and Subsequent Short Term Optimization for the Single Owner Model.
    We are interested in how the single owner model performs compared to the water regulation Procedure.
    Therefore, we need to extract the water used, the hourly bidding curves, etc.
    What is more efficient, who earns the most money?
    
"""
function SingleOwnerSimulation(R::Vector{Reservoir}, K::Vector{HydropowerPlant}, mu_up, mu_down, inflow_data, price_data, MediumModelSingle, currentweek::Int64, price_point_count::Int64, iteration_count_bidding::Int64, T::Int64, stage_count_bidding::Int64, stage_count_short::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64)
    _, f = AverageReservoirLevel(R, inflow_data)
    cuts =  ReservoirLevelCutsSingle(R, K,  f, currentweek, stage_count_short) 
    WaterCuts = WaterValueCutsSingle(R, MediumModelSingle, cuts, currentweek)
    Ω, P, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_bidding)
    PPoints = Test_Price_Points(Ω, scenario_count_prices, T, mu_up)
    HourlyBiddingCurve = SingleOwnerBidding(R, K, PPoints, Ω, P , cuts, WaterCuts, mu_up, mu_down, iteration_count_bidding, T, stage_count_bidding)
    # price = Price_Scenarios_Short(price_data, 1, stage_count_short)[1][1]
    price = Ω[stage_count_bidding][rand(1:scenario_count_prices)].price
    inflow = Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, scenario_count_inflows)[1]
    Obligations = MarketClearingSolo(price,  HourlyBiddingCurve, PPoints, T)
    Qnom, z_up, z_down = SingleOwnerScheduling(R, K, Obligations, price, inflow, Ω, P, cuts, WaterCuts, mu_up, mu_down, iteration_count_short, T, stage_count_short)
    println("current reservoir: $([(r.dischargepoint, r.currentvolume) for r in R])")
    update_reservoir!(Qnom, Dict(r => inflow[r][1] for r in R))
    println("after operation: $([(r.dischargepoint, r.currentvolume) for r in R])")
    return HourlyBiddingCurve, Obligations, price, Qnom, z_up, z_down
end

"""
Simulate The Entire Procedure for one week, and obtain the solution afterwards. We are interested in
    - Total revenue
    - Hourly Bidding Bidding Curves
    - Day Ahead Obligations
    - Realized Price
    - Reservoir Level Changes (Real / Individual)
    """
function ExampleSimulation(R::Vector{Reservoir}, J::Vector{Participant}, mu_up, mu_down, inflow_data, price_data, MediumModel_j, MediumModel_O, currentweek::Int64, scenario_count_prices::Int64; ColumnReservoir = ColumnReservoir)
    # PPoints = Create_Price_Points(price_data, price_point_count)
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)

    cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
    Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
    cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
    WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
    WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)

    HourlyBiddingCurves, Qnoms1, Ω1, PPoints = FirstLayerSimulation(J, R, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, iteration_count_short, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices)
    # price = Price_Scenarios_Short(price_data, 1, stage_count_short)[1][1]
    price = Ω1[J[1]][stage_count_bidding][rand(1:scenario_count_prices)].price
    inflow = Dict(r => Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, scenario_count_inflows)[1][r][1] for r in R)
    println(inflow)
    Obligations = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
    Qadj1, _, P_Swap1, _, _, _ = water_regulation(Qnoms1, Qref, inflow, false)

    Qnoms2 = SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligations, price, price_data, inflow_data, Qref, cuts,  WaterCuts, iteration_count_short, mu_up, mu_down, T, stage_count_short)
    Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)

    z_ups, z_downs = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligations, mu_up, mu_down, T)
    return HourlyBiddingCurves, Obligations, Qnoms1, Qadj1, P_Swap1, price, Qnoms2, Qadj2, P_Swap2, z_ups, z_downs
end


"""
    FirstLayerSimulation(Strategy)

    We simulate the Bidding Models which are relevant depending on the participants strategy.
    We include the generation of the nonanticipatory/anticipatory uncertainty sets here
    The HourlyBiddingCurves and Nominations are returned.
"""
function FirstLayerSimulation(J::Vector{Participant}, all_res::Vector{Reservoir}, Strategy::Dict{Participant, String}, price_data::DataFrame, inflow_data::DataFrame, Qref::Dict{Reservoir, Float64}, cuts, cutsOther, WaterCuts, WaterCutsOther, iteration_count_short::Int64, mu_up::Float64, mu_down::Float64, T::Int64, stage_count_bidding::Int64, scenario_count_prices::Int64)

    HourlyBiddingCurves = Dict{Participant, Dict{Int64, Vector{Float64}}}()
    Qnoms = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    Ω1 = Dict{Participant, Any}()
    PPoints = Dict{Participant, Vector{Vector{Float64}}}()
    for j in J
        R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
        Ω_NA_local, P_NA_local, Ω_scenario_local, P_scenario_local = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, all_res, stage_count_bidding)
        PPoints[j] = Test_Price_Points(Ω_NA_local, scenario_count_prices, T, mu_up)
        println(PPoints)
        if Strategy[j] == "Nonanticipatory"
            Qnom, HourlyBiddingCurve = Nonanticipatory_Bidding(R, j, PPoints[j], Ω_NA_local, P_NA_local, Qref, cuts[j], WaterCuts[j], iteration_count_short, mu_up, mu_down, T, stage_count_bidding)
            HourlyBiddingCurves[j] = HourlyBiddingCurve
        else
            Ω_A_local, P_A_local = create_Ω_Anticipatory(Ω_NA_local, Ω_scenario_local, P_scenario_local, J, j, all_res, PPoints[j], Qref, cutsOther, WaterCutsOther, stage_count_bidding, mu_up, mu_down, T)
            println(Ω_A_local)
            Qnom, HourlyBiddingCurve = Anticipatory_Bidding(all_res, j, J, PPoints[j], Ω_A_local, P_A_local, Qref, cuts[j], WaterCuts[j], iteration_count_short, mu_up, mu_down, T, stage_count_bidding)
            HourlyBiddingCurves[j] = HourlyBiddingCurve
        end
        for r in all_res
            if j.participationrate[r] > 0
                Qnoms[(participant = j, reservoir = r)] = Qnom[r]
            else
                Qnoms[(participant = j, reservoir = r)] = Qref[r]
            end
        end
        Ω1[j] = Ω_NA_local
    end
    return HourlyBiddingCurves, Qnoms, Ω1, PPoints
end

"""

    SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligations, price_data, inflow_data, Qref, cuts, WaterCuts, iteration_count_short, mu_up, mu_down, T, stage_count_short)

    Simulate the 
"""
function SecondLayerSimulation(J::Vector{Participant}, R::Vector{Reservoir}, Qnoms1, Qadj1, Obligations, price, price_data::DataFrame, inflow_data::DataFrame, Qref, cuts, WaterCuts, iteration_count_short::Int64, mu_up::Float64, mu_down::Float64, T::Int64, stage_count_short::Int64)
    
    QnomO1 = OthersNomination(Qnoms1, Qadj1, J, R)
    Qnoms2 = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    for j in J
        Ω_NA_local, P_NA_local, _, _ = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_short)
        Qnom = ShortTermScheduling(R, j, J, Qref, Obligations[j], price, QnomO1[j], Ω_NA_local, P_NA_local, cuts[j], WaterCuts[j], iteration_count_short, mu_up,mu_down, T, stage_count_short) 
        for r in R
            if j.participationrate[r] > 0
                Qnoms2[(participant = j, reservoir = r)] = Qnom[r]
            else
                Qnoms2[(participant = j, reservoir = r)] = Qref[r]
            end
        end
    end

    return Qnoms2
end

"""
ThirdLayerSimulation()

    Final Layer of Decision making for all Participants.
    With fixed Power Swaps and Obligations, make a Decision about how to schedule the power plants over the day.
    Necessary to find realized revenue for all participants.
"""

function ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligations, mu_up, mu_down, T)
    z_ups = Dict{Participant, Vector{Float64}}()
    z_downs = Dict{Participant, Vector{Float64}}()
    for j in J
        _, z_up, z_down, w = RealTimeBalancing(
            R,
            j,
            Qadj2,
            P_Swap2[j],
            T,
            mu_up,
            mu_down,
            Obligations[j])
        z_ups[j] = z_up
        z_downs[j] = z_down
    end
    return z_ups, z_downs
end

# HourlyBiddingCurves, Obligations, Qnoms1, Qadj1, P_Swap1, price, Qnoms2, Qadj2, P_Swap2, z_ups, z_downs = ExampleSimulation(R, J, mu_up, mu_down,inflow_data, price_data, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded, currentweek, scenario_count_prices)

HourlyBiddingCurve, Obligation, price_solo, Qnom, z_up, z_down = SingleOwnerSimulation(R, K, mu_up, mu_down, inflow_data, price_data, MediumModelSingle, currentweek, price_point_count, iteration_count_Bidding, T, stage_count_bidding, stage_count_short, scenario_count_prices, scenario_count_inflows)



# Final_Revenue(J, price, Obligations, z_ups, z_downs, mu_up, mu_down, T)
Final_Revenue_Solo(price_solo, Obligation, z_up, z_down, mu_up, mu_down, T)