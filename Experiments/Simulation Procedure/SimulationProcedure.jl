
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
const savepath_results = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\SingleVsIndividual\\SingleVsIndividual.csv"
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
    SingleOwnerSimulation(R, K, price_data, inflow_data, cuts, WaterCuts, PPoints, iteration_count, mu_up, mu_down, T, stage_count_bidding)

    Simulate the Bidding and Subsequent Short Term Optimization for the Single Owner Model.
    We are interested in how the single owner model performs compared to the water regulation Procedure.
    Therefore, we need to extract the water used, the hourly bidding curves, etc.
    What is more efficient, who earns the most money?
    
"""
function SingleOwnerSimulation(R::Vector{Reservoir}, K::Vector{HydropowerPlant}, mu_up, mu_down, inflow_data, price_data, MediumModelSingle, currentweek::Int64, iteration_count_bidding::Int64, T::Int64, stage_count_bidding::Int64, stage_count_short::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64)
    _, f = AverageReservoirLevel(R, inflow_data)
    cuts =  ReservoirLevelCutsSingle(R, K,  f, currentweek, stage_count_short) 
    WaterCuts = WaterValueCutsSingle(R, MediumModelSingle, cuts, currentweek)
    Ω, P, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_bidding)
    PPoints = Create_Price_Points(Ω, scenario_count_prices, T, mu_up)
    HourlyBiddingCurve = SingleOwnerBidding(R, K, PPoints, Ω, P , cuts, WaterCuts, mu_up, mu_down, iteration_count_bidding, T, stage_count_bidding)
    # price = Price_Scenarios_Short(price_data, 1, stage_count_short)[1][1]
    price = Ω[stage_count_bidding][rand(1:scenario_count_prices)].price
    inflow = Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, scenario_count_inflows)[1]
    Obligations = MarketClearingSolo(price,  HourlyBiddingCurve, PPoints, T)
    Qnom, z_up, z_down = SingleOwnerScheduling(R, K, Obligations, price, inflow, Ω, P, cuts, WaterCuts, mu_up, mu_down, iteration_count_short, T, stage_count_short)
    # update_reservoir!(Qnom, Dict(r => inflow[r][1] for r in R))
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
    Obligations = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
    Qadj1, _, P_Swap1, _, _, _ = water_regulation(Qnoms1, Qref, inflow, false)

    Qnoms2 = SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligations, price, price_data, inflow_data, Qref, cuts,  WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_short)
    Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)

    z_ups, z_downs = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligations, mu_up, mu_down, T)
    return HourlyBiddingCurves, Obligations, Qnoms1, Qadj1, P_Swap1, price, Qnoms2, Qadj2, P_Swap2, z_ups, z_downs
end

"""
function Comparison_Single_VF(R, J, mu_up, mu_down, inflow_data, price_data, MediumModel_j, MediumModel_O, currentweek, scenario_count_prices)

    Experiment to compare the performance of the Single Owner Model and the Producers within the Water Regulation Procedure
    Every model more or less works independently, though the realized prices and inflow affect both systems equally: We compare the simulations over a week under equal realizations
    We are interested in the reservoir levels, revenues, and especially the ratio between both.
    Furthermore, we would like to know how the individual reservoirs behave, and if there are participants that profit differently from other players in the system. How would their profits have been,
    if they were distributed from the single owner model based on total production power?
"""
function Comparison_Single_VF(R::Vector{Reservoir},J::Vector{Participant}, mu_up::Float64, mu_down::Float64, inflow_data::DataFrame, price_data::DataFrame,
    MediumModel_j, MediumModel_O, MediumModelSingle, Initial_Reservoir::Dict{Reservoir, Float64}, Initial_Individual_Reservoir::Dict{Participant, Dict{Reservoir, Float64}}, currentweek::Int64, scenario_count_prices::Int64, scenario_count_inflows::Int64, iteration_count_bidding::Int64, Strategy::Dict{Participant, String}; ColumnReservoir = ColumnReservoir)
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Individual_Revenues = Dict{Participant, Float64}[]
    Revenues_single = Float64[]
    l_singles = Dict{Reservoir, Float64}[]
    l_vfs = Dict{Reservoir, Float64}[]
    l_vf_inds = Dict{Participant, Dict{Reservoir, Float64}}[]
    Qadjs = Dict{Reservoir, Float64}[]
    Qnoms = Dict{Reservoir, Float64}[]
    Qnoms_individual = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}[]
    P_Swaps = Dict{Participant, Dict{Reservoir, Float64}}[]
    Obligations = Dict{Participant, Vector{Float64}}[]
    Obligations_solo = Vector{Float64}[]
    z_ups = Dict{Participant, Float64}[]
    z_downs = Dict{Participant, Float64}[]
    l_single = Dict(r => Initial_Reservoir[r] for r in R)
    for r in R
        r.currentvolume = Initial_Reservoir[r]
        for j in J
            j.individualreservoir[r] = Initial_Individual_Reservoir[j][r]
        end
    end
    for day in 1:days
        println("current day: ", day)
        # Reservoir Parameters
        Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
        # Water Value Function(s)
        cuts_single =  ReservoirLevelCutsSingle(R, K,  f, currentweek, stage_count_short) 
        WaterCuts_single = WaterValueCutsSingle(R, MediumModelSingle, cuts_single, currentweek)
        cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
        Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
        cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
        WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModel_j[j], cuts[j], currentweek) for j in J)
        WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModel_O[j], cutsOther[j], currentweek) for j in J)
        # Bidding Problem: Single Owner and Indiviual Problems
        Ω_single, P_single, _, _= create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_bidding)
        PPoints_single = Create_Price_Points(Ω_single, scenario_count_prices, T, mu_up)
        HourlyBiddingCurve = SingleOwnerBidding(R, l_single, K, PPoints_single, Ω_single, P_single , cuts_single, WaterCuts_single, mu_up, mu_down, iteration_count_bidding, T, stage_count_bidding)
        HourlyBiddingCurves, Qnoms1, Ω1, PPoints = FirstLayerSimulation(J, R, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek)
        # Realization of uncertain parameters: inflow and price
        price = create_Ω_Nonanticipatory(price_data, inflow_data, 1, 1, currentweek, R, stage_count_bidding)[1][stage_count_bidding][1].price
        inflow_array = Inflow_Scenarios_Short(inflow_data, currentweek, R, stage_count_short, 1)[1]
        inflow = Dict(r => inflow_array[r][1] for r in R)
        # Market Clearing
        Obligation_solo = MarketClearingSolo(price, HourlyBiddingCurve, PPoints_single, T)
        Obligation = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
        # Single Owner Scheduling
        Qnom, z_up_solo, z_down_solo = SingleOwnerScheduling(R, l_single, K, Obligation_solo, price, inflow_array, Ω_single, P_single, cuts_single, WaterCuts_single, mu_up, mu_down, iteration_count_short, T, stage_count_short)
        
        # Individual Scheduling: Two Optimization problems and water regulation procedure (includes reservoir update)
        Qadj1, _, P_Swap1, _, _, _ = water_regulation(Qnoms1, Qref, inflow, false)
        Qnoms2 = SecondLayerSimulation(J, R, Qnoms1, Qadj1, Obligation, price, price_data, inflow_data, Qref, cuts,  WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_short, scenario_count_prices, scenario_count_inflows, currentweek)
        Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, inflow, true)
        z_up, z_down = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligation, mu_up, mu_down, T)
        #Reservoir Update (Single Owner)
        Initial_Reservoir = Dict{Reservoir, Float64}(r => r.currentvolume for r in R)
        for r in R
            l_single[r] = l_single[r] + inflow[r] - Qnom[r]
            for j in J
                Initial_Individual_Reservoir[j][r] = j.individualreservoir[r] - Qnoms2[(participant = j, reservoir = r)] + Qref[r]
            end
        end
        # Revenues
        Individual_Revenue = Final_Revenue(J, price, Obligation, z_up, z_down, mu_up, mu_down, T)
        Revenue_single = Final_Revenue_Solo(price, Obligation_solo, z_up_solo, z_down_solo, mu_up, mu_down, T)
        println(l_single)
        push!(Individual_Revenues, Individual_Revenue)
        push!(Revenues_single, Revenue_single)
        push!(Qnoms, Qnom)
        push!(Qnoms_individual, Qnoms2)
        push!(Obligations, Obligation)
        push!(Obligations_solo, Obligation_solo)
        push!(P_Swaps, P_Swap2)
        push!(Qadjs, Qadj2)
        push!(l_singles, Dict(r => l_single[r] for r in R))
        push!(l_vfs, Initial_Reservoir)
        push!(l_vf_inds, Dict(j => Dict(r =>  Initial_Individual_Reservoir[j][r] for r in R) for j in J))
        push!(z_ups, Dict(j => sum(z_up[j][t] for t in 1:T) for j in J))
        push!(z_downs, Dict(j => sum(z_down[j][t] for t in 1:T) for j in J))
    end  
    return Individual_Revenues, Revenues_single, l_singles, l_vfs, l_vf_inds, Qadjs, Qnoms, Qnoms_individual, P_Swaps, Obligations, Obligations_solo, z_ups, z_downs
end

# HourlyBiddingCurves, Obligations, Qnoms1, Qadj1, P_Swap1, price, Qnoms2, Qadj2, P_Swap2, z_ups, z_downs = ExampleSimulation(R, J, mu_up, mu_down,inflow_data, price_data, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded, currentweek, scenario_count_prices)
# HourlyBiddingCurve, Obligation, price_solo, Qnom, z_up, z_down = SingleOwnerSimulation(R, K, mu_up, mu_down, inflow_data, price_data, MediumModelSingle, currentweek, iteration_count_Bidding, T, stage_count_bidding, stage_count_short, scenario_count_prices, scenario_count_inflows)

"""
function SplitRevenues(J, R, Revenues)

    Based on the revenues generated from the single owner model, split up the revenues according to the participationrate

"""
function SplitRevenues(J::Vector{Participant}, R::Vector{Reservoir}, Revenues::Vector{Float64})::Vector{Dict{Participant, Float64}}
    TotalParticipationRates = Dict{Participant, Float64}(j => sum(j.participationrate[r] for r in R)/(sum(p.participationrate[r] for r in R for p in J)) for j in J)
    println(TotalParticipationRates)
    SplitRevenue = [Dict{Participant, Float64}(j => TotalParticipationRates[j] * Revenues[t] for j in J) for t in 1:2]
    return SplitRevenue
end

# Final_Revenue(J, price, Obligations, z_ups, z_downs, mu_up, mu_down, T)
# Final_Revenue_Solo(price_solo, Obligation, z_up, z_down, mu_up, mu_down, T)

Strategy = Dict(j => "Nonanticipatory" for j in J)

Weeks = [15]
WeeklyAverageReservoirLevels = Dict(week => Dict(r => mean(AverageReservoirLevel(R, inflow_data)[1][r][(week-1)*7 + 1: week * 7]) for r in R) for week in 1:52)
Individual_Revenues = Dict{Int64, Vector{Dict{Participant, Float64}}}()
Revenues_single = Dict{Int64, Vector{Float64}}()
l_singles = Dict{Int64, Vector{Dict{Reservoir, Float64}}}()
l_vfs = Dict{Int64, Vector{Dict{Reservoir, Float64}}}()
l_vf_inds = Dict{Int64, Vector{Dict{Participant, Dict{Reservoir, Float64}}}}()
Qadjs = Dict{Int64, Vector{Dict{Reservoir, Float64}}}()
Qnoms = Dict{Int64, Vector{Dict{Reservoir, Float64}}}()
Qnoms_individual = Dict{Int64, Vector{Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}}()
P_Swaps = Dict{Int64, Vector{Dict{Participant, Dict{Reservoir, Float64}}}}()
Split_Revenues = Dict{Int64, Vector{Dict{Participant, Float64}}}()
Obligations = Dict{Int64, Vector{Dict{Participant, Vector{Float64}}}}()
Obligations_solo = Dict{Int64, Vector{Vector{Float64}}}()
z_ups = Dict{Int64, Vector{Dict{Participant, Float64}}}()
z_downs = Dict{Int64, Vector{Dict{Participant, Float64}}}()
for week in Weeks
    println("Current week: $(week)")
    Initial_Reservoir = copy(WeeklyAverageReservoirLevels[week])
    Initial_Individual_Reservoir = Dict{Participant, Dict{Reservoir, Float64}}(j => copy(WeeklyAverageReservoirLevels[week]) for j in J)
    Individual_Revenue, Revenue_single, l_single, l_vf, l_vf_ind, Qadj, Qnom, Qnom_individual, P_Swap, Obligation, Obligation_solo, z_up, z_down = Comparison_Single_VF(R, J, mu_up, mu_down, inflow_data, price_data, MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded, MediumModelSingle, Initial_Reservoir, Initial_Individual_Reservoir, week, scenario_count_prices, scenario_count_inflows, iteration_count_bidding, Strategy)
    Split_Revenue = SplitRevenues(J, R, Revenue_single)
    Individual_Revenues[week] = Individual_Revenue
    Revenues_single[week] = Revenue_single
    l_singles[week] = l_single
    l_vfs[week] = l_vf
    l_vf_inds[week] = l_vf_ind
    Qadjs[week] = Qadj
    Qnoms[week] = Qnom
    Qnoms_individual[week] = Qnom_individual
    P_Swaps[week] = P_Swap
    Split_Revenues[week] = Split_Revenue
    Obligations[week] = Obligation
    Obligations_solo[week] = Obligation_solo
    z_ups[week] = z_up
    z_downs[week] = z_down
end
"""
function ResultsToDataFrameSingleVsIndividual()
    
    To save the results for later analysis, organize them inside a DataFrame.
    This is also helpful to do some statisical analysis, with functions from DataFrames.jl
"""
function ResultsToDataFrameSingleVsIndividual(savepath, J, R, Strategy, Individual_Revenues, Revenue_single, l_single, l_vf, l_vf_ind, Qadj, Qsingle, Qnoms_individual, P_Swaps, Obligations, Obligations_single, z_ups, z_downs, Weeks::Vector{Int64}; save = true)
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


df_results = ResultsToDataFrameSingleVsIndividual(savepath_results, J, R, Strategy, Individual_Revenues, Revenues_single, l_singles, l_vfs, l_vf_inds, Qadjs, Qnoms, Qnoms_individual, P_Swaps, Obligations, Obligations_solo, z_ups, z_downs, Weeks)