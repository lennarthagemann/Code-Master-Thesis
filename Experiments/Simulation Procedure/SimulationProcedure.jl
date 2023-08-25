
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
using Serialization
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

const ColumnReservoir = Dict{Reservoir, String}(R[1] => "Flasjon Inflow", R[2] => "Holmsjon Inflow")
const scenario_count_inflows = 1
const scenario_count_prices = 3
const stage_count_short = 3
const stage_count_bidding = 2
const price_point_count = 5
const T = 24
const currentweek = 2
const iteration_count_short = 50
const iteration_count_Bidding = 100

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)
NameToParticipant = Dict{String, Participant}(j.name => j for j in J)
MediumModelDictionary_j_loaded = ReadMediumModel(savepath_watervalue * "\\Participant.jls", NameToParticipant);
MediumModelDictionary_O_loaded = ReadMediumModel(savepath_watervalue * "\\OtherParticipant.jls", NameToParticipant);

Strategy = Dict(j => "Anticipatory" for j in J)


"""
Simulate The Entire Procedure for one week, and obtain the solution afterwards. We are interested in
    - Total revenue
    - Hourly Bidding Bidding Curves
    - Day Ahead Obligations
    - Realized Price
    - Reservoir Level Changes (Real / Individual)
"""
function ExampleSimulation(R::Vector{Reservoir}, J::Vector{Participant}, inflow_data, price_data, MediumModel_j, MediumModel_O, currentweek::Int64, price_point_count::Int64; ColumnReservoir = ColumnReservoir)
    PPoints = Create_Price_Points(price_data, price_point_count)
    mu_up, mu_down = BalanceParameters(price_data)
    
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)

    
    cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
    Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
    cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
    WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
    WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)

    HourlyBiddingCurves, Qnoms1 = FirstLayerSimulation(J, R, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, PPoints, iteration_count_short, mu_up, mu_down, T, stage_count_bidding, stage_count_short)

    return HourlyBiddingCurves, Qnoms1
end
"""
    SingleOwnerSimulation(R, K, price_data, inflow_data, cuts, WaterCuts, PPoints, iteration_count, mu_up, mu_down, T, stage_count_bidding)

    Simulate the Bidding and Subsequent Short Term Optimization for the Single Owner Model.
    We are interested in how the single owner model performs compared to the water regulation Procedure.
    Therefore, we need to extract the water used, the hourly bidding curves, etc.
    What is more efficient, who earns the most money?
    
"""
function SingleOwnerSimulation()
    
end

"""
    FirstLayerSimulation(Strategy)

    We simulate the Bidding Models which are relevant depending on the participants strategy.
    We include the generation of the nonanticipatory/anticipatory uncertainty sets here
    The HourlyBiddingCurves and Nominations are returned.
"""
function FirstLayerSimulation(J::Vector{Participant}, all_res::Vector{Reservoir}, Strategy::Dict{Participant, String}, price_data::DataFrame, inflow_data::DataFrame, Qref::Dict{Reservoir, Float64}, cuts, cutsOther, WaterCuts, WaterCutsOther, PPoints::Vector{Float64}, iteration_count_short::Int64, mu_up::Float64, mu_down::Float64, T::Int64, stage_count_short::Int64, stage_count_bidding::Int64)

    HourlyBiddingCurves = Dict{Participant, Dict{Int64, Vector{Float64}}}()
    Qnoms = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()

    for j in J
        R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
        Ω_NA_local, P_NA_local, Ω_scenario_local, P_scenario_local = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, all_res, stage_count_bidding, ColumnReservoir)
        if Strategy[j] == "Nonanticipatory"
            println(cuts[j])
            Qnom, HourlyBiddingCurve = Nonanticipatory_Bidding(R, j, PPoints, Ω_NA_local, P_NA_local, Qref, cuts[j], WaterCuts[j], iteration_count_short, mu_up, mu_down, T, stage_count_bidding)
            HourlyBiddingCurves[j] = HourlyBiddingCurve
        else
            println(R)
            Ω_A_local, P_A_local = create_Ω_Anticipatory(Ω_NA_local, Ω_scenario_local, P_scenario_local, J, j, all_res, PPoints, Qref, cutsOther, WaterCutsOther, stage_count_bidding, mu_up, mu_down, T)
            Qnom, HourlyBiddingCurve = Anticipatory_Bidding(all_res, j, J, PPoints, Ω_A_local, P_A_local, Qref, cuts[j], WaterCuts[j], iteration_count_short, mu_up, mu_down, T, stage_count_bidding)
            HourlyBiddingCurves[j] = HourlyBiddingCurve
        end
        for r in all_res
            if j.participationrate[r] > 0
                Qnoms[(participant = j, reservoir = r)] = Qnom[r]
            else
                Qnoms[(participant = j, reservoir = r)] = Qref[r]
            end
        end
    end
    return HourlyBiddingCurves, Qnoms
end

HourlyBiddingCurves, Qnoms = ExampleSimulation(R, J, inflow_data, price_data, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded, currentweek, price_point_count)