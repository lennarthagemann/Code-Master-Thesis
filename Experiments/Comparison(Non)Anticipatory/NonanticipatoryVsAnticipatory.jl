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
using Tables
import Base
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
const savepath_experiment ="C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\NonanticipatoryVsAnticipatory"
R, K, J = read_data(filepath_Ljungan)

const ColumnReservoir = Dict(r => r.dischargepoint * " Inflow" for r in R)
const scenario_count_inflows = 1
const scenario_count_prices = 20
const scenario_count_prices_medium = 3
const scenario_count_inflows_weekly = 3
const stage_count_short = 2
const stage_count_bidding = 2
const stage_count_medium = 52
const price_point_count = 5
const T = 24
const iteration_count_short = 20
const iteration_count_bidding = 10
const iteration_count_medium = 1000

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
    Initial_Reservoir::Dict{Reservoir, Float64}, Initial_Individual_Reservoir::Dict{Participant, Dict{Reservoir, Float64}}, MediumModel_j::Dict{Participant, SDDP.PolicyGraph{Int64}}, MediumModel_O::Dict{Participant, SDDP.PolicyGraph{Int64}},
    currentweek::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, iteration_count_bidding::Int64, iteration_count_short::Int64)
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
    cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
    cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
    WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
    WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)
    Qnoms_Bidding = Dict{Dict{Participant, String}, Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}()
    Qnoms_Scheduling = Dict{Dict{Participant, String}, Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}()
    Qadjs = Dict{Dict{Participant, String}, Dict{Reservoir, Float64}}()
    P_Swaps = Dict{Dict{Participant, String}, Dict{Participant, Dict{Reservoir, Float64}}}()
    z_ups = Dict{Dict{Participant, String}, Dict{Participant, Vector{Float64}} }()
    z_downs = Dict{Dict{Participant, String}, Dict{Participant, Vector{Float64}}}()
    Individual_Revenues = Dict{Dict{Participant, String}, Dict{Participant, Float64}}()
    l_reals = Dict{Dict{Participant, String}, Dict{Reservoir, Float64}}()
    l_inds = Dict{Dict{Participant, String}, Dict{Participant, Dict{Reservoir, Float64}}}()
    Obligations = Dict{Dict{Participant, String}, Dict{Participant, Vector{Float64}}}()
    for strat in Strategy_Combinations
        println("Current Strat: \n$(strat)")
        for r in R
            r.currentvolume = Initial_Reservoir[r]
            for j in J
                j.individualreservoir[r] = Initial_Individual_Reservoir[j][r]
            end
        end
        HourlyBiddingCurves, Qnoms1, Ω1, PPoints = FirstLayerSimulation(J, R, strat, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_bidding, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek; printlevel = 0)
        price = Ω1[J[1]][stage_count_bidding][rand(1:scenario_count_prices)].price
        inflow = Dict(r => Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, scenario_count_inflows)[1][r][1] for r in R)
        Obligation = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
        Qadj1, _, P_Swap1, _, _, _ = water_regulation(Qnoms1, Qref, inflow, false)
        Qnoms2 = SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligation, price, price_data, inflow_data, Qref, cuts,  WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_short, scenario_count_prices, scenario_count_inflows, currentweek)
        Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)
        z_up, z_down = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligation, mu_up, mu_down, T)
        Individual_Revenue = Final_Revenue(J, price, Obligation, z_up, z_down, mu_up, mu_down, T)
        
        # for j in J
        #     for r in R
        #         l_vf_ind[j][r] = l_vf_ind[j][r] - Qnoms2[(participant = j, reservoir = r)] + Qref[r]
        #     end
        # end
        Obligations[strat] = Obligation
        Qnoms_Bidding[strat] = Qnoms1
        Qnoms_Scheduling[strat] = Qnoms2
        Qadjs[strat] = Qadj2
        P_Swaps[strat] = P_Swap2
        z_ups[strat] = z_up
        z_downs[strat] = z_down
        Individual_Revenues[strat] = Individual_Revenue
        l_reals[strat] = Dict{Reservoir, Float64}(r => r.currentvolume for r in R)
        l_inds[strat] = Dict{Participant, Dict{Reservoir, Float64}}(j => Dict(r => j.individualreservoir[r] for r in R) for j in J)
    end

    return Qnoms_Bidding, Obligations, Qnoms_Scheduling, Qadjs, P_Swaps, z_ups, z_downs, Individual_Revenues, l_reals, l_inds
end


"""
function ResultsToDataFrame()
    
    To save the results for later analysis, organize them inside a DataFrame.
    This is also helpful to do some statisical analysis, with functions from DataFrames.jl
    """
function ResultsToDataFrame(savepath, J::Vector{Participant}, R::Vector{Reservoir}, Strategy_Combinations, Qnoms_Bidding, Obligations, Qnoms_Scheduling, Qadjs, P_Swaps, z_ups, z_downs, Individual_Revenues, l_reals, l_inds, currentweek::Int64; save = true)
    column_names_df_nominations = ["week", ["Strategy_" * j.name for j in J]..., ["Qnom1_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qnom2_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., ["P_Swap_" * j.name * "_" * r.dischargepoint for j in J for r in R]...]
    column_types_df_nominations = [Int64, [String for j in J]..., [Float64 for j in J for r in R]..., [Float64 for j in J for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]...]
    column_names_df_Obligations = ["week", ["Strategy_" * j.name for j in J]..., ["Obligations_" * j.name for j in J]..., ["z_up_" * j.name for j in J]..., ["z_down_" *j.name for j in J]..., ["Revenue_"*j.name for j in J]...]
    column_types_df_Obligations = [Int64, [String for j in J]..., [Vector{Float64} for j in J]..., [Vector{Float64} for j in J]..., [Vector{Float64} for j in J]..., [Float64 for j in J]...]
    column_names_df_Reservoirs = ["week", ["Strategy_" * j.name for j in J]...,  ["l_real_" * r.dischargepoint for r in R]..., ["l_ind_" * j.name * "_" * r.dischargepoint for j in J for r in R]...]
    column_types_df_Reservoirs = [Int64, [String for j in J]...,  [Float64 for r in R]..., [Float64 for j in J for r in R]...,]
    
    if isfile(savepath * "\\NominationsBounded.csv")
        # File exists, attempt to load the DataFrame from the file
        df_nominations = CSV.File(savepath * "\\NominationsBounded.csv", types = column_types_df_nominations) |> DataFrame
        df_Obligations = CSV.File(savepath * "\\ObligationsBounded.csv", types = column_types_df_Obligations) |> DataFrame
        df_Reservoirs = CSV.File(savepath * "\\ReservoirsBounded.csv", types = column_types_df_Reservoirs) |> DataFrame
        println("DataFrame already exists. Add input parameters as new data...")
        println(eltype.(eachcol(df_nominations)))
        @assert names(df_nominations) == column_names_df_nominations
        @assert (eltype.(eachcol(df_nominations))) == column_types_df_nominations
        @assert names(df_Obligations) == column_names_df_Obligations
        @assert (eltype.(eachcol(df_Obligations))) == column_types_df_Obligations
        @assert names(df_Reservoirs) == column_names_df_Reservoirs
        @assert (eltype.(eachcol(df_Reservoirs))) == column_types_df_Reservoirs
    else
        # File doesn't exist
        println("DataFrame does not exist yet and will subsequently be created and filled...")
        df_nominations = DataFrame()
        df_Obligations = DataFrame()
        df_Reservoirs = DataFrame()
        for (name, type) in zip(column_names_df_nominations, column_types_df_nominations)
            df_nominations[!, name] = Vector{type}()
        end
        for (name, type) in zip(column_names_df_Obligations, column_types_df_Obligations)
            df_Obligations[!, name] = Vector{type}()
        end
        for (name, type) in zip(column_names_df_Reservoirs, column_types_df_Reservoirs)
            df_Reservoirs[!, name] = Vector{type}()
        end
        println(names(df_nominations))
        println(names(df_Obligations))
        println(names(df_Reservoirs))
    end

    for strat in Strategy_Combinations
        println(strat)
        row_nominations = (week = currentweek,
        Strategy_Sydkraft = strat[J[1]],
        Strategy_Fortum = strat[J[2]],
        Strategy_Statkraft = strat[J[3]],
        Qnom1_Sydkraft_Flasjon = Qnoms_Bidding[strat][(participant = J[1], reservoir = R[1])],
        Qnom1_Sydkraft_Holmsjon = Qnoms_Bidding[strat][(participant = J[1], reservoir = R[2])],
        Qnom1_Fortum_Flasjon = Qnoms_Bidding[strat][(participant = J[2], reservoir = R[1])],
        Qnom1_Fortum_Holmsjon = Qnoms_Bidding[strat][(participant = J[2], reservoir = R[2])],
        Qnom1_Statkraft_Flasjon = Qnoms_Bidding[strat][(participant = J[3], reservoir = R[1])],
        Qnom1_Statkraft_Holmsjon = Qnoms_Bidding[strat][(participant = J[3], reservoir = R[2])],
        Qnom2_Sydkraft_Flasjon = Qnoms_Scheduling[strat][(participant = J[1], reservoir = R[1])],
        Qnom2_Sydkraft_Holmsjon = Qnoms_Scheduling[strat][(participant = J[1], reservoir = R[2])],
        Qnom2_Fortum_Flasjon = Qnoms_Scheduling[strat][(participant = J[2], reservoir = R[1])],
        Qnom2_Fortum_Holmsjon = Qnoms_Scheduling[strat][(participant = J[2], reservoir = R[2])],
        Qnom2_Statkraft_Flasjon = Qnoms_Scheduling[strat][(participant = J[3], reservoir = R[1])],
        Qnom2_Statkraft_Holmsjon = Qnoms_Scheduling[strat][(participant = J[3], reservoir = R[2])],
        Qadj_Flasjon = Qadjs[strat][R[1]],
        Qadj_Holmsjon = Qadjs[strat][R[2]],
        P_Swap_Sydkraft_Flasjon = P_Swaps[strat][J[1]][R[1]],
        P_Swap_Sydkraft_Holmsjon = P_Swaps[strat][J[1]][R[2]],
        P_Swap_Fortum_Flasjon = P_Swaps[strat][J[2]][R[1]],
        P_Swap_Fortum_Holmsjon = P_Swaps[strat][J[2]][R[2]],
        P_Swap_Statkraft_Flasjon = P_Swaps[strat][J[3]][R[1]],
        P_Swap_Statkraft_Holmsjon = P_Swaps[strat][J[3]][R[2]]
        )
        row_Obligations = (week = currentweek,
        Strategy_Sydkraft = strat[J[1]],
        Strategy_Fortum = strat[J[2]],
        Strategy_Statkraft = strat[J[3]],
        Obligations_Sydkraft = Obligations[strat][J[1]],
        Obligations_Fortum = Obligations[strat][J[2]],
        Obligations_Statkraft = Obligations[strat][J[3]],
        z_up_Sydkraft = z_ups[strat][J[1]],
        z_up_Fortum = z_ups[strat][J[2]],
        z_up_Statkraft = z_ups[strat][J[3]],
        z_down_Sydkraft = z_downs[strat][J[1]],
        z_down_Fortum = z_downs[strat][J[2]],
        z_down_Statkraft = z_downs[strat][J[3]],
        Revenue_Sydkraft = Individual_Revenues[strat][J[1]],
        Revenue_Fortum = Individual_Revenues[strat][J[2]],
        Revenue_Statkraft = Individual_Revenues[strat][J[3]]
        )
        row_Reservoirs = (week = currentweek,
        Strategy_Sydkraft = strat[J[1]],
        Strategy_Fortum = strat[J[2]],
        Strategy_Statkraft = strat[J[3]],
        l_real_Flasjon = l_reals[strat][R[1]],
        l_real_Holmsjon = l_reals[strat][R[2]],
        l_ind_Sydkraft_Flasjon = l_inds[strat][J[1]][R[1]],
        l_ind_Sydkraft_Holmsjon = l_inds[strat][J[1]][R[2]],
        l_ind_Fortum_Flasjon = l_inds[strat][J[2]][R[1]],
        l_ind_Fortum_Holmsjon = l_inds[strat][J[2]][R[2]],
        l_ind_Statkraft_Flasjon = l_inds[strat][J[3]][R[1]],
        l_ind_Statkraft_Holmsjon = l_inds[strat][J[3]][R[2]],
        )
        push!(df_nominations, row_nominations)
        push!(df_Obligations, row_Obligations)
        push!(df_Reservoirs, row_Reservoirs)
    end

    if save == true
        println("Results will be saved at $(savepath)...")
        println(df_Obligations)
        CSV.write(savepath * "\\NominationsBounded.csv", df_nominations)
        CSV.write(savepath * "\\ObligationsBounded.csv", df_Obligations)
        CSV.write(savepath * "\\ReservoirsBounded.csv", df_Reservoirs)
    end
    return df_nominations, df_Obligations, df_Reservoirs
end

Weeks = [20, 30, 40, 45, 50]
WeeklyAverageReservoirLevels = Dict(week => Dict(r => mean(AverageReservoirLevel(R, inflow_data)[1][r][(week-1)*7 + 1: week * 7]) for r in R) for week in 1:52)
for i in 1:50
    for week in Weeks
        currentweek = week
        Initial_Reservoir = WeeklyAverageReservoirLevels[currentweek]
        Initial_Individual_Reservoir = Dict{Participant, Dict{Reservoir, Float64}}(j => WeeklyAverageReservoirLevels[currentweek] for j in J)
        Qnoms_Bidding, Obligations, Qnoms_Scheduling, Qadjs, P_Swaps, z_ups, z_downs, Individual_Revenues, l_reals, l_inds = AnticipatoryVsNonanticipatory(R, J, mu_up, mu_down, inflow_data, price_data,
        Initial_Reservoir, Initial_Individual_Reservoir, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded, currentweek, scenario_count_prices, scenario_count_inflows, iteration_count_bidding, iteration_count_short)
        df_nominations, df_Obligations, df_Reservoirs = ResultsToDataFrame(savepath_experiment, J, R, Strategy_Combinations, Qnoms_Bidding, Obligations, Qnoms_Scheduling, Qadjs, P_Swaps, z_ups, z_downs, Individual_Revenues, l_reals, l_inds, currentweek)
    end
end