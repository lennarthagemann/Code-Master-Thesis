
#=
---------------------------------------------------------------------------------------------

It is interesting to see, how much profit the producers would have acquired, if they were alone
in the river. 
How does it differ in comparison to the situation in the water regulation setup?
For that we simulate both models for the producer (single and VF) and compare the results.
 
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


# Overload the parse function to handle Vector{Float64}
Base.parse(::Type{Vector{Float64}}, s::AbstractString) = parse_vector_string(s)
Base.tryparse(::Type{Vector{Float64}}, s::AbstractString) = parse_vector_string(s)

# Custom function to parse the string into a Vector{Float64}
function parse_vector_string(s::AbstractString)
    values_str = replace(s, r"[\[\]]" => "")
    values = split(values_str, ",")

    # Convert the string elements to Float64 and create a Vector{Float64}
    vector_obj = parse.(Float64, values)
    return vector_obj
end


const filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
const filepath_prices = pwd() * "\\Inflow Forecasting\\Data\\Spot Prices\\prices_df.csv"
const filepath_inflows = pwd() * "\\Inflow Forecasting\\Data\\Inflow\\Data from Flasjoen and Holmsjoen.csv"
const savepath_watervalue = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Inflow Forecasting\\WaterValue"
const savepath_results = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\SingleVsIndividual\\SingleVsIndividualBoundedNom.csv"
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
const days = 7
const currentweek = 50
const iteration_count_short = 50
const iteration_count_bidding = 10
const iteration_count_medium = 1000

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)

PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)
Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded = ReadMediumModel(savepath_watervalue, J, R, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
MediumModelSingle = ReadMediumModelSingle(savepath_watervalue, R, K, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
mu_up, mu_down = BalanceParameters(price_data)

Strategy = Dict(j => "Nonanticipatory" for j in J)
"""

function ComparisonRealSingleOwner()

    Comparison of Waster Regulation Procedure to Single Owner Revenues side by side.
    Done simultaneuously for all participants.

"""

function ComparisonRealSingleOwner(J::Vector{Participant},
    all_res::Vector{Reservoir},
    price_data::DataFrame, inflow_data::DataFrame,
    WeeklyAverageReservoirLevels::Dict{Int64, Dict{Reservoir, Float64}}, MediumModel_j, MediumModel_O,
    mu_up::Float64, mu_down::Float64, T::Int64, stage_count_bidding::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, currentweek::Int64; iteration_count = iteration_count_bidding, printlevel = 0, stability_report = false, bounded_bidding = true)
    
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    HourlyBiddingCurvesVF = Dict{Participant, Dict{Int64, Vector{Float64}}}[]
    HourlyBiddingCurvesSingle = Dict{Participant, Dict{Int64, Vector{Float64}}}[]
    QnomsVF = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}[]
    Qnoms_single = Dict{Participant, Dict{Reservoir, Float64}}[]
    z_ups = Dict{Participant, Vector{Float64}}[]
    z_downs = Dict{Participant, Vector{Float64}}[]
    z_ups_solo = Dict{Participant, Vector{Float64} }[]
    z_downs_solo = Dict{Participant, Vector{Float64}}[]
    Individual_RevenuesVF = Dict{Participant, Float64}[]
    Individual_Revenues_solo = Dict{Participant, Float64}[]
    Initial_Reservoir = Dict(r => WeeklyAverageReservoirLevels[currentweek][r] for r in R)
    l_single = Dict(j => Dict(r => Initial_Reservoir[r] for r in R) for j in J)
    Initial_Individual_Reservoir = Dict{Participant, Dict{Reservoir, Float64}}(j => copy(WeeklyAverageReservoirLevels[currentweek]) for j in J)
    for r in all_res
        r.currentvolume = Initial_Reservoir[r]
    end
    for day in 1:days
        println("current day: ", day)
        # Reservoir Parameters
        Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
        # Water Value Functions for all Participants
        cuts = Dict(j => ReservoirLevelCuts(all_res, j.plants, j, f, currentweek, stage_count_short) for j in J)
        WaterCuts = Dict(j => WaterValueCuts(all_res, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
        Others = Dict(j => OtherParticipant(J,j,all_res)[1] for j in J)
        cutsOther = Dict(j => ReservoirLevelCuts(all_res, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
        WaterCutsOther = Dict(j => WaterValueCuts(all_res, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)
        # Bidding for all Participants
        println(typeof(Initial_Reservoir))
        HourlyBiddingCurveVF, QnomVF, _, PPointsVF = FirstLayerSimulation(J, K, all_res, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_bidding, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek)
        # Bidding for all Participants in simulated solo environment
        HourlyBiddingCurve, PPoints = SingleOwnerParticipantBidding(J, all_res, price_data, inflow_data, l_single, cuts, WaterCuts, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek)
        # Market Clearing
        price = create_Ω_Nonanticipatory(price_data, inflow_data, 1, 1, currentweek, all_res, stage_count_bidding)[1][stage_count_bidding][1].price
        inflow_array = Inflow_Scenarios_Short(inflow_data, currentweek, all_res, stage_count_short, 1)[1]
        inflow = Dict(r => inflow_array[r][1] for r in all_res)
        ObligationVF = MarketClearing(price, HourlyBiddingCurveVF, PPointsVF, J, T)
        Obligation = MarketClearing(price, HourlyBiddingCurve, PPoints, J, T)
        # Scheduling for all Participants
        # Single Owner Scheduling
        Ω_single, P_single, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, all_res, stage_count_bidding)
        Qnom, z_up_solo, z_down_solo = SingleOwnerParticipantScheduling(J, all_res, Initial_Reservoir, Obligation, price, inflow_array, cuts, WaterCuts, mu_up, mu_down, T, stage_count_short, scenario_count_prices, scenario_count_inflows, currentweek)
        # Water Regulation and Single Owner Scheduling VF
        Qadj1, _, _ , _, _, _ = water_regulation(QnomVF, Qref, inflow, false)
        Qnoms2 = SecondLayerSimulation(J, R, QnomVF, Qadj1, ObligationVF, price, price_data, inflow_data, Qref, cuts,  WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_short, scenario_count_prices, scenario_count_inflows, currentweek)
        # Real Time Balancing for all Participants VF
        Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)
        z_up, z_down = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, ObligationVF, mu_up, mu_down, T)
        #Reservoir Update (Single Owner)
        Initial_Reservoir =  Dict{Reservoir, Float64}(r => r.currentvolume for r in R)
        for j in J
            for r in collect(filter(r -> j.participationrate[r] > 0, all_res))
                l_single[j][r] = l_single[j][r] + inflow[r] - Qnom[j][r]
                Initial_Individual_Reservoir[j][r] = j.individualreservoir[r] - Qnoms2[(participant = j, reservoir = r)] + Qref[r]
            end
        end
        # Final Revenues
        Individual_RevenueVF = Final_Revenue(J, price, ObligationVF, z_up, z_down, mu_up, mu_down, T)
        Individual_Revenue_solo = Final_Revenue(J, price, Obligation, z_up, z_down, mu_up, mu_down, T)
        # Save in Dictionary
        push!(HourlyBiddingCurvesVF, HourlyBiddingCurveVF)
        push!(HourlyBiddingCurvesSingle, HourlyBiddingCurve)
        push!(z_ups, z_up)
        push!(z_downs, z_down)
        push!(z_ups_solo, z_up_solo)
        push!(z_downs_solo, z_down_solo)
        push!(Individual_RevenuesVF, Individual_RevenueVF)
        push!(Individual_Revenues_solo, Individual_Revenue_solo)
        push!(QnomsVF, QnomVF)
        push!(Qnoms_single, Qnom)
    end
    return HourlyBiddingCurvesVF, HourlyBiddingCurvesSingle, QnomsVF, Qnoms_single, z_ups, z_downs, z_ups_solo, z_downs_solo, Individual_RevenuesVF, Individual_Revenues_solo
end

"""

function SingleOwnerParticipantBidding()
    
    Do an equivalent Single Owner evaluation of the bidding problem.
    To be used alongside the Bidding within Water Regulation context.
    Returns the Bidding Curves for all participants

"""
function SingleOwnerParticipantBidding(J::Vector{Participant},
    all_res::Vector{Reservoir},
    price_data::DataFrame, inflow_data::DataFrame,
    l_single::Dict{Participant, Dict{Reservoir, Float64}}, cuts, WaterCuts,
    mu_up::Float64, mu_down::Float64, T::Int64, stage_count_bidding::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, currentweek::Int64; iteration_count = iteration_count_bidding, printlevel = 0, stability_report = false, bounded_bidding = true)
    
    _ , f = AverageReservoirLevel(R, inflow_data)
    HourlyBiddingCurves = Dict{Participant, Dict{Int64, Vector{Float64}}}()
    for j in J
        println("Current Participant: $(j.name)")
        R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
        println("Relevant Reservoirs: $(R)")
        Ω_single, P_single, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_bidding)
        PPoints_single = Create_Price_Points(Ω_single, scenario_count_prices, T, mu_up)
        HourlyBiddingCurves[j] = SingleOwnerBidding(R, l_single[j], j.plants, PPoints_single, Ω_single, P_single, cuts[j], WaterCuts[j], mu_up, mu_down, iteration_count, T, stage_count_bidding)
        PPoints[j] = PPoints_single
    end
    for j in J
        for t in 1:T
            sort!(HourlyBiddingCurves[j][t])
        end
    end
    return HourlyBiddingCurves, PPoints
end

"""

function SingleOwnerParticipantScheduling()
    
    Do an equivalent Single Owner evaluation of the bidding problem.
    To be used alongside the Bidding within Water Regulation context

"""
function SingleOwnerParticipantScheduling(J::Vector{Participant},
    all_res::Vector{Reservoir}, Initial_Reservoir,
    y_initial, price, f_initial,
    cuts, WaterCuts,
    mu_up::Float64, mu_down::Float64, T::Int64, stage_count_short::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, currentweek::Int64; iteration_count = iteration_count_short)
    Qnoms = Dict{Participant, Dict{Reservoir, Float64}}()
    z_ups = Dict{Participant, Vector{Float64}}()
    z_downs = Dict{Participant, Vector{Float64}}()
    for j in J
        R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
        Ω_local, P_local, _, _ = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_short)
        Qnom, z_up, z_down = SingleOwnerScheduling(R, Initial_Reservoir, j.plants,
        y_initial[j], price, f_initial,
        Ω_local, P_local, cuts[j], WaterCuts[j], mu_up, mu_down, iteration_count, T, stage_count_short)
        Qnoms[j] = Qnom
        z_ups[j] = z_up
        z_downs[j] = z_down
    end 
    return Qnoms, z_ups, z_downs
end

"""

function ResultsToDataFrameRealSingleOwner()
    
    To save the results for later analysis, organize them inside a DataFrame.
    This is also helpful to do some statisical analysis, with functions from DataFrames.jl

"""
function ResultsToDataFrameRealSingleOwner(savepath, J, R, Strategy, QnomsVF, Qnoms_single, Individual_RevenuesVF, Individual_Revenues_solo, Weeks::Vector{Int64}, days::Int64; save = true)
    column_names = ["week", "day", ["Strategy_" * j.name for j in J]..., ["Qnoms_VF_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qnoms_single_" * j.name * "_" * r.dischargepoint for  j in J for r in R]..., ["Individual_Revenues_VF_" * j.name for j in J]..., ["Individual_Revenues_single_" * j.name for j in J]...]
    column_types = [Int64, Int64, [String for j in J]..., [Float64 for j in J for r in R]..., [Float64 for j in J for r in R]..., [Float64 for j in J]..., [Float64 for j in J]...]
    @assert length(column_names) == length(column_types)
    if isfile(savepath)
        df = CSV.File(savepath, types = column_types) |> DataFrame 
    else
        df = DataFrame()
        for (name, type) in zip(column_names, column_types)
            df[!, name] = Vector{type}()
        end
    end
    for week in Weeks
        for day in 1:days
            row = (
                week = week,
                day = day,
                Strategy_Sydkraft = Strategy[J[1]],
                Strategy_Fortum = Strategy[J[2]],
                Strategy_Statkraft = Strategy[J[3]],
                Qnoms_single_Sydkraft_Flasjon = Qnoms_single[week][day][J[1]][R[1]],
                Qnoms_single_Sydkraft_Holmsjon = Qnoms_single[week][day][J[1]][R[2]],
                Qnoms_single_Fortum_Flasjon = Qnoms_single[week][day][J[2]][R[1]],
                Qnoms_single_Fortum_Holmsjon = Qnoms_single[week][day][J[2]][R[2]],
                Qnoms_single_Statkraft_Flasjon = Qnoms_single[week][day][J[3]][R[1]],
                Qnoms_single_Statkraft_Holmsjon = Qnoms_single[week][day][J[3]][R[2]],
                Qnoms_VF_Sydkraft_Flasjon = QnomsVF[week][day][(participant = J[1], reservoir = R[1])], 
                Qnoms_VF_Sydkraft_Holmsjon = QnomsVF[week][day][(participant = J[1], reservoir = R[2])], 
                Qnoms_VF_Fortum_Flasjon = QnomsVF[week][day][(participant = J[2], reservoir = R[1])], 
                Qnoms_VF_Fortum_Holmsjon = QnomsVF[week][day][(participant = J[2], reservoir = R[2])], 
                Qnoms_VF_Statkraft_Flasjon = QnomsVF[week][day][(participant = J[3], reservoir = R[1])], 
                Qnoms_VF_Statkraft_Holmsjon = QnomsVF[week][day][(participant = J[3], reservoir = R[2])], 
                VF_Revenues_Sydkraft = Individual_RevenuesVF[week][day][J[1]],
                VF_Revenues_Fortum = Individual_RevenuesVF[week][day][J[2]],
                VF_Revenues_Statkraft = Individual_RevenuesVF[week][day][J[3]],
                Single_Revenue = Revenue_single[week][day],
                Split_Revenues_Sydkraft = Individual_Revenues_solo[week][day][J[1]],
                Split_Revenues_Fortum = Individual_Revenues_solo[week][day][J[2]],
                Split_Revenues_Statkraft = Individual_Revenues_solo[week][day][J[3]],
            )
            println(length(names(df)))
            push!(df, row)
        end
    end
    if save == true
        println("DataFrame will be saved at $(savepath)...")
        CSV.write(savepath, df)
    end
    return df
end


Weeks = [2]
# Weeks = [1, 14, 19, 29, 39, 44, 49]
# Weeks = [9, 13, 18, 28, 38, 43, 48]
# Weeks = [8, 12, 17, 27, 37, 42, 47]
# Weeks = [7, 11, 16, 26, 36, 41, 46]


# _, f = AverageReservoirLevel(R, inflow_data)
# Initial_Reservoir = WeeklyAverageReservoirLevels[2]
# HourlyBiddingCurvesSingle, PPoints = SingleOwnerParticipantBidding(J, R, price_data, inflow_data, WeeklyAverageReservoirLevels[2], MediumModelDictionary_j_loaded, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, 2)

# price = create_Ω_Nonanticipatory(price_data, inflow_data, 1, 1, currentweek, R, stage_count_bidding)[1][stage_count_bidding][1].price
# inflow_array = Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, 1)[1]
# inflow = Dict(r => inflow_array[r][1] for r in R)
# Obligation = MarketClearing(price, HourlyBiddingCurvesSingle, PPoints, J, T)

# cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
# WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModelDictionary_j_loaded[j], cuts[j], currentweek) for j in J)

# Qnoms, z_ups, z_downs = SingleOwnerParticipantScheduling(J, R, Initial_Reservoir, Obligation, price, inflow_array, cuts, WaterCuts, mu_up, mu_down, T, stage_count_short, scenario_count_prices, scenario_count_inflows, currentweek;)


WeeklyAverageReservoirLevels = Dict(week => Dict(r => mean(AverageReservoirLevel(R, inflow_data)[1][r][(week-1)*7 + 1: week * 7]) for r in R) for week in 1:52)
QnomsVF = Dict{Int64, Vector{Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}}() 
Qnoms_single = Dict{Int64, Vector{Dict{Participant, Dict{Reservoir, Float64}}}}()
z_ups = Dict{Int64, Vector{Dict{Participant, Vector{Float64}}}}() 
z_downs = Dict{Int64, Vector{Dict{Participant, Vector{Float64}}}}() 
z_ups_solo = Dict{Int64, Vector{Dict{Participant, Vector{Float64}}}}() 
z_downs_solo = Dict{Int64, Vector{Dict{Participant, Vector{Float64}}}}() 
Individual_RevenuesVF = Dict{Int64, Vector{Dict{Participant, Float64}}}() 
Individual_Revenues_solo = Dict{Int64, Vector{Dict{Participant, Float64}}}() 
for week in Weeks
    println("----------- current week: $(week) --------------")
    _, _, QnomVF, Qnom_single, z_up, z_down, z_up_solo, z_down_solo, Individual_RevenueVF, Individual_Revenue_solo = ComparisonRealSingleOwner(J, R, price_data, inflow_data,
        WeeklyAverageReservoirLevels, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded,
        mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek;)
    QnomsVF[week] = QnomVF
    Qnoms_single[week] = Qnom_single
    z_ups[week] = z_up
    z_downs[week] = z_down
    z_ups_solo[week] = z_up_solo
    z_downs_solo[week] = z_down_solo
    Individual_RevenuesVF[week] = Individual_RevenueVF
    Individual_Revenues_solo[week] = Individual_Revenue_solo
end
        


df_results = ResultsToDataFrameSingleVsIndividual(savepath_results, J, R, Strategy, QnomsVF, Qnoms_single, Individual_RevenuesVF, Individual_Revenues_solo, Weeks)