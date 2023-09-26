
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
const iteration_count_bidding = 100
const iteration_count_medium = 1000

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)

PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)
Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded = ReadMediumModel(savepath_watervalue, J, R, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
MediumModelSingle = ReadMediumModelSingle(savepath_watervalue, R, K, Ω_medium, P_medium, stage_count_medium, iteration_count_medium)
mu_up, mu_down = BalanceParameters(price_data)

"""
function ComparisonRealSingleOwner()

    Comparison of Waster Regulation Procedure to Single Owner Revenues side by side.
    Done simultaneuously for all participants.


"""
function ComparisonRealSingleOwner(J::Vector{Participant})
    return
end

"""
function SingleOwnerParticipantBidding()
    
    Do an equivalent Single Owner evaluation of the bidding problem.
    To be used alongside the Bidding within Water Regulation context

"""
function SingleOwnerParticipantBidding(J::Vector{Participant},
    all_res::Vector{Reservoir},
    price_data::DataFrame, inflow_data::DataFrame,
    cuts,  WaterCuts,
    Initial_Reservoir::Dict{Reservoir, Float64}, 
    mu_up::Float64, mu_down::Float64, T::Int64, stage_count_bidding::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, currentweek::Int64; printlevel = 0, stability_report = false, bounded_bidding = true)

    HourlyBiddingCurves = Dict{Participant, Dict{Int64, Vector{Float64}}}()
    for j in J
        R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
        Ω_single, P_single, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_bidding)
        PPoints_single = Create_Price_Points(Ω_single, scenario_count_prices, T, mu_up)
        HourlyBiddingCurves[j] = SingleOwnerBidding(R, Initial_Reservoir, j.plants, PPoints_single, Ω_single, P_single, cuts[j], WaterCuts[j], mu_up, mu_down, iteration_count, T)
    end
    return HourlyBiddingCurves
end

"""
function SingleOwnerParticipantScheduling()
    
    Do an equivalent Single Owner evaluation of the bidding problem.
    To be used alongside the Bidding within Water Regulation context

"""
function SingleOwnerParticipantScheduling(J::Vector{Participant})
    for j in J
        Qnom, z_up, z_down = SingleOwnerScheduling(R, Initial_Reservoir, j.plants,
        y_initial, price, f_initial,
        Ω, P, cuts, WaterCuts, mu_up, mu_down,iteration_count, T, Stages)
        Qnoms[j] = Qnom
        z_ups[j] = z_up
        z_downs[j] = z_down
    end 
end


"""
function ResultsToDataFrameRealSingleOwner()
    
    To save the results for later analysis, organize them inside a DataFrame.
    This is also helpful to do some statisical analysis, with functions from DataFrames.jl
"""
function ResultsToDataFrameRealSingleOwner(savepath, J, R, Strategy, Individual_Revenues, Revenue_single, l_single, l_vf, l_vf_ind, Qadj, Qsingle, Qnoms_individual, P_Swaps, Obligations, Obligations_single, z_ups, z_downs, Weeks::Vector{Int64}; save = true)
    column_names = ["week", "day", ["Strategy_" * j.name for j in J]..., ["Reservoir_Single_" * r.dischargepoint for r in R]..., ["Reservoir_VF_" * r.dischargepoint for r in R]..., ["Individual_Reservoir_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., ["Q_single_" * r.dischargepoint for r in R]..., ["Qnom_$(j.name)_$(r.dischargepoint)" for j in J for r in R]..., ["P_Swaps_$(j.name)_$(r.dischargepoint)" for j in J for r in R]..., ["Obligations_" * j.name for j in J]...,"Obligations_single", ["VF_Revenues_" * j.name for j in J]..., "Single_Revenue", ["Split_Revenues_" * j.name for j in J]..., ["z_ups_" * j.name for j in J]..., ["z_downs_" * j.name for j in J]... ]
    column_types = [Int64, Int64, [String for j in J]..., [Float64 for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]..., [Float64 for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]..., [Float64 for j in J for r in R]..., [Vector{Float64} for j in J]..., Vector{Float64},  [Float64 for j in J]..., Float64, [Float64 for j in J]..., [Float64 for j in J]..., [Float64 for j in J]...]
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
                Reservoir_Single_Flasjon = l_single[week][day][R[1]],
                Reservoir_Single_Holmsjon = l_single[week][day][R[2]],
                Reservoir_VF_Flasjon = l_vf[week][day][R[1]],
                Reservoir_VF_Holmsjon = l_vf[week][day][R[2]],
                Individual_Reservoir_Sydkraft_Flasjon = l_vf_ind[week][day][J[1]][R[1]],
                Individual_Reservoir_Sydkraft_Holmsjon = l_vf_ind[week][day][J[1]][R[2]],
                Individual_Reservoir_Fortum_Flasjon = l_vf_ind[week][day][J[2]][R[1]],
                Individual_Reservoir_Fortum_Holmsjon = l_vf_ind[week][day][J[2]][R[2]],
                Individual_Reservoir_Statkraft_Flasjon = l_vf_ind[week][day][J[3]][R[1]],
                Individual_Reservoir_Statkraft_Holmsjon = l_vf_ind[week][day][J[3]][R[2]],
                Qadj_Flasjon = Qadj[week][day][R[1]],
                Qadj_Holmsjon = Qadj[week][day][R[2]],
                Q_single_Flasjon = Qsingle[week][day][R[1]],
                Q_single_Holmsjon = Qsingle[week][day][R[2]],
                Qnom_Sydkraft_Flasjon = Qnoms_individual[week][day][(participant = J[1], reservoir = R[1])], 
                Qnom_Sydkraft_Holmsjon = Qnoms_individual[week][day][(participant = J[1], reservoir = R[2])], 
                Qnom_Fortum_Flasjon = Qnoms_individual[week][day][(participant = J[2], reservoir = R[1])], 
                Qnom_Fortum_Holmsjon = Qnoms_individual[week][day][(participant = J[2], reservoir = R[2])], 
                Qnom_Statkraft_Flasjon = Qnoms_individual[week][day][(participant = J[3], reservoir = R[1])], 
                Qnom_Statkraft_Holmsjon = Qnoms_individual[week][day][(participant = J[3], reservoir = R[2])], 
                P_Swaps_Sydkraft_Flasjon = P_Swaps[week][day][J[1]][R[1]], 
                P_Swaps_Sydkraft_Holmsjon = P_Swaps[week][day][J[1]][R[2]], 
                P_Swaps_Fortum_Flasjon = P_Swaps[week][day][J[2]][R[1]], 
                P_Swaps_Fortum_Holmsjon = P_Swaps[week][day][J[2]][R[2]], 
                P_Swaps_Statkraft_Flasjon = P_Swaps[week][day][J[3]][R[1]], 
                P_Swaps_Statkraft_Holmsjon = P_Swaps[week][day][J[3]][R[2]], 
                Obligations_Sydkraft = Obligations[week][day][J[1]],
                Obligations_Fortum = Obligations[week][day][J[2]],
                Obligations_Statkraft = Obligations[week][day][J[3]],
                Obligations_single = Obligations_single[week][day],
                VF_Revenues_Sydkraft = Individual_Revenues[week][day][J[1]],
                VF_Revenues_Fortum = Individual_Revenues[week][day][J[2]],
                VF_Revenues_Statkraft = Individual_Revenues[week][day][J[3]],
                Single_Revenue = Revenue_single[week][day],
                Split_Revenues_Sydkraft = Split_Revenues[week][day][J[1]],
                Split_Revenues_Fortum = Split_Revenues[week][day][J[2]],
                Split_Revenues_Statkraft = Split_Revenues[week][day][J[3]],
                z_ups_Sydkraft = z_ups[week][day][J[1]],
                z_ups_Fortum = z_ups[week][day][J[2]],
                z_ups_Statkraft = z_ups[week][day][J[3]],
                z_downs_Sydkraft = z_downs[week][day][J[1]],
                z_downs_Fortum = z_downs[week][day][J[2]],
                z_downs_Statkraft = z_downs[week][day][J[3]]
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


# df_results = ResultsToDataFrameSingleVsIndividual(savepath_results, J, R, Strategy, Individual_Revenues, Revenues_single, l_singles, l_vfs, l_vf_inds, Qadjs, Qnoms, Qnoms_individual, P_Swaps, Obligations, Obligations_solo, z_ups, z_downs, Weeks)