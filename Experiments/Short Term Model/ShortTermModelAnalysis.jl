#=
-----------------------------------------------------------------------------------------------

The ShortTermModel is used to find the final nomination in the Water Regulation Procedure.
It is based on the other participants' nomination, as well as the realized price and obligation.
We are interested on what parameters the model is sensitive. Most importantly:

    - How sensitive is the model towards Others' nomination, if price and inflow is fixed?

It should not interfere too significantly with my own nomination, as a huge power swap will already be
dealt with through the additional water at disposal. Some complications are expected though, especially if the adjusted flow causes spillage.
Then the problem becomes balancing between spillage, power swap and obligations to fulfil.
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
const savepath_results = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results"

R, K, J = read_data(filepath_Ljungan)

const ColumnReservoir = Dict(r => r.dischargepoint * " Inflow" for r in R)
const scenario_count_inflows = 1
const scenario_count_prices = 10
const scenario_count_prices_medium = 3
const scenario_count_inflows_weekly = 3
const stage_count_short = 2
const stage_count_bidding = 2
const stage_count_medium = 52
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
mu_up, mu_down = BalanceParameters(price_data)
Strategy = Dict{Participant, String}(j => "Nonanticipatory" for j in J)

"""
function create_Obligation( T, J, R, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, iteration_count_short, mu_up, mu_down, stage_count_bidding, scenario_count_prices)

    Helper function to create Obligation to be used in ShortTermProblem.
"""
function create_Obligations(T::Int64, J::Vector{Participant}, R::Vector{Reservoir}, Strategy::Dict{Participant, String}, price_data::DataFrame, inflow_data::DataFrame, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir::Dict{Reservoir, Float64}, Initial_Individual_Reservoir::Dict{Participant, Dict{Reservoir, Float64}}, iteration_count_short, mu_up, mu_down, stage_count_bidding, scenario_count_prices, scenario_count_inflows, current_week)
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    HourlyBiddingCurves, Qnoms, Ω, PPoints = FirstLayerSimulation(J, R, Strategy, price_data, inflow_data, Qref, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, T, stage_count_bidding, scenario_count_prices, scenario_count_inflows, current_week)
    price = create_Ω_Nonanticipatory(price_data, inflow_data, 1, 1, currentweek, R, stage_count_bidding)[1][stage_count_bidding][1].price
    Obligation = MarketClearing(price, HourlyBiddingCurves, PPoints, J, T)
    return Obligation, price, Qnoms
end
"""
function create_Qadjs(R, K, count)

    Helper Function to decide what kind of adjusted flow are to be used in the ShortTerScheduling Analysis.
    count gives the amount of different nominations to run through the short term scheduling
"""
function create_Qadjs(R::Vector{Reservoir}, K::Vector{HydropowerPlant}, count::Int64)::Vector{Dict{Reservoir, Float64}}
    Qadjs = Dict{Reservoir, Float64}[]
    max_flow = max([k.spillreference for k in K]...)
    vals = [c * 1/count * max_flow for c in 1:count]
    for val in vals
        for c in 1:count
            Qadj = Dict{Reservoir, Float64}(R[1] => val * (c/count), R[2] => val * (1 - c/count))
            push!(Qadjs, Qadj)
        end
    end
    return Qadjs
end

"""
function create_Qadjs_simple(R, K, count)

    Helper Function to decide what kind of adjusted flow are to be used in the ShortTerScheduling Analysis.
    count gives the amount of different nominations to run through the short term scheduling.
    In comparison to the other function, the adjusted flow only varies around the nomination given.
"""
function create_Qadjs_simple(R::Vector{Reservoir}, QnomsDict::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, inflow_data::DataFrame, count::Int64, currentweek::Int64)::Dict{Participant, Vector{Dict{Reservoir, Float64}}}
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    QadjsDict = Dict{Participant, Vector{Dict{Reservoir, Float64}}}()
    vals = Dict{Participant, Vector{Float64}}()
    for j in J
        Qadjs = Dict{Reservoir, Float64}[]
        if j.name == "Sydkraft"
            MeanTotalFlow = QnomsDict[(participant = j, reservoir = R[1])]
            DeviationTotalFlow  = MeanTotalFlow/4
            val = collect(range(MeanTotalFlow - DeviationTotalFlow, MeanTotalFlow + DeviationTotalFlow, count^2))
            vals[j] = val
        else
            MeanTotalFlow = sum(QnomsDict[(participant = j, reservoir = r)] for r in R)
            DeviationTotalFlow  = MeanTotalFlow/4
            val = collect(range(MeanTotalFlow - DeviationTotalFlow, MeanTotalFlow + DeviationTotalFlow, count))
            vals[j] = val
        end
        if j.name == "Sydkraft"
            for val in vals[j]
                Qadj = Dict{Reservoir, Float64}(R[1] => val, R[2] => Qref[R[2]])
                push!(Qadjs, Qadj)
            end
        else
            for val in vals[j]
                for c in 1:count
                    Qadj = Dict{Reservoir, Float64}(R[1] => (c-1)/count * val, R[2] =>( 1 - (c-1)/count) * val)
                    push!(Qadjs, Qadj)
                end
            end
        end
        QadjsDict[j] = Qadjs
    end
    return QadjsDict
end
"""
function create_simple_Obligation(J)

    Helper function to create Obligation to be used in ShortTermProblem.
    In this case we look at an Obligation, that is relatively smooth and easy to handle, even if a power swap has to be delivered
    => Generate an Obligation, that is at about 50% of the maximum spill reference of the power plants owned, plus a little bit of variation (+- 20%)

"""
function create_simple_Obligations(J::Vector{Participant}, R::Vector{Reservoir}, T::Int64, price_data::DataFrame, inflow_data::DataFrame, currentweek::Int64, stage_count_bidding::Int64)
    price = create_Ω_Nonanticipatory(price_data, inflow_data, 1, 1, currentweek, R, stage_count_bidding)[1][stage_count_bidding][1].price
    Obligations = Dict{Participant, Vector{Float64}}()
    Qnoms = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    for j in J
        K = j.plants
        AvgObligation = sum(k.equivalent for k in K) * 0.5 * min([k.spillreference for k in K]...)
        Deviation = sum(k.equivalent for k in K) * 0.2 * min([k.spillreference for k in K]...)
        dist = Uniform(-Deviation, Deviation)
        Obligations[j] = [AvgObligation + rand(dist) for t in 1:T]
        Obligations[j][1] = 0
        for r in R
            if j.participationrate[r] > 0
                # Only Sydkraft nominates in one reservoir, the rest nominate in two reservoirs.
                if j.name == "Sydkraft"
                    Qnoms[(participant = j, reservoir = r)] = sum(Obligations[j][t] for t in 1:T)/(T * j.participationrate[r])
                else
                    Qnoms[(participant = j, reservoir = r)] = sum(Obligations[j][t] for t in 1:T)/(T * j.participationrate[r] * 2)
                end
            else
                Qnoms[(participant = j, reservoir = r)] = Qref[r]
            end
        end
    end
    return Obligations, price, Qnoms
end

Obligations_simple, price_simple, Qnoms_simple = create_simple_Obligations(J, R, T, price_data, inflow_data, currentweek, stage_count_bidding)

"""
function ShortTermProblemParametrizedNomination(R, K, j, Obligation, price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, stage_count_short)

    The point of this function is to execute the ShortTermProblem for different input parameters. Mainly:

        -Different nomination by other people while Day Ahead Obligations are the same
    
    This leads to different power swap situations -> Technically this should be the same, but having to deliver a power swap might in some cases lead 
    to higher/lower nominations. The constraint of having to deliver a power swap at certaint times together with obligations
    might make it impossible to fulfill an obligation at a certain time.

    The Obligations will be taken from real obligations, resulting from either low/middle/high Bidding Curves.

    It will also be interesting to see if the participants are interested in the reservoirs differently, even if the water volumes
    reaching their power plants are the same.


"""
function ShortTermProblemParametrized(R::Vector{Reservoir}, K::Vector{HydropowerPlant},J::Vector{Participant}, j::Participant, Qnom_initial::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Obligations::Dict{Participant, Vector{Float64}}, cuts, WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, mu_up, mu_down,  price_data::DataFrame, inflow_data::DataFrame, scenario_count_prices::Int64, scenario_count_inflows::Int64, currentweek::Int64, stage_count_short::Int64, T::Int64, count_adjusted_flows::Int64)
    println("Current Participant: $(j.name)")
    Qadjs = create_Qadjs(R, K, count_adjusted_flows)
    Qnoms =  Dict{Dict{Reservoir, Float64}, Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}()
    QnomOs = Dict{Dict{Reservoir, Float64}, Dict{Participant, Dict{Reservoir, Float64}}}()
    QadjsFinal = Dict{Dict{Reservoir, Float64}, Dict{Reservoir, Float64}}()
    z_ups = Dict{Dict{Reservoir, Float64}, Float64}()
    z_downs = Dict{Dict{Reservoir, Float64}, Float64}()
    Revenues = Dict{Dict{Reservoir, Float64}, Float64}()
    l_traj, f = AverageReservoirLevel(R, inflow_data)
    Qref = CalculateReferenceFlow(R, l_traj, f, currentweek)
    for Qadj in Qadjs
        println("Current Adjusted Flow before Scheduling: $(Qadj)")
        QnomO = OthersNomination(Qnom_initial, Qadj, J, R)
        Qnoms2 = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
        Ω_NA_local, P_NA_local, _, _ = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_short)
        Qnom = ShortTermScheduling(R, j, J, Qref, Obligations[j], price, QnomO[j], Ω_NA_local, P_NA_local, cuts[j], WaterCuts[j], iteration_count_short, mu_up,mu_down, T, stage_count_short, Initial_Reservoir, Initial_Individual_Reservoir; printlevel = 0) 
        for r in R
            if j.participationrate[r] > 0
                Qnoms2[(participant = j, reservoir = r)] = Qnom[r]
            else
                Qnoms2[(participant = j, reservoir = r)] = Qref[r]
            end
            for notj in filter(p -> p.name != j.name,J)
                Qnoms2[(participant = notj, reservoir = r)] = QnomO[j][r]
            end
        end
        Qadj2, _, P_Swap2, _, _, _ = water_regulation(Qnoms2, Qref, Dict(r => f[r][currentweek] for r in R), false)
        println("Nomination after Scheduling: $(Qnom)")
        println("Adjusted Flow after Scheduling: $(Qadj2)")
        z_up, z_down = ThirdLayerSimulation(J, R, Qadj2, P_Swap2, Obligations, mu_up, mu_down, T)
        Revenue  = Final_Revenue(J, price, Obligations, z_up, z_down, mu_up, mu_down, T)
        Qnoms[Qadj] = Qnoms2
        QnomOs[Qadj] = QnomO
        QadjsFinal[Qadj] = Qadj2
        z_ups[Qadj] = sum(z_up[j][t] for t in 1:T)
        z_downs[Qadj] =  sum(z_down[j][t] for t in 1:T)
        Revenues[Qadj] = Revenue[j]
    end
    return Qnoms, QnomOs, QadjsFinal, z_ups, z_downs, Revenues
end


l_traj, f = AverageReservoirLevel(R, inflow_data)
WeeklyAverageReservoirLevels = Dict(week => Dict(r => mean(AverageReservoirLevel(R, inflow_data)[1][r][(week-1)*7 + 1: week * 7]) for r in R) for week in 1:52)
cuts = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, currentweek, stage_count_short) for j in J)
Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, currentweek, stage_count_short) for j in J)
WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModelDictionary_j_loaded[j], cuts[j], currentweek) for j in J)
WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModelDictionary_O_loaded[j], cutsOther[j], currentweek) for j in J)
Initial_Reservoir = WeeklyAverageReservoirLevels[currentweek]
Initial_Individual_Reservoir = Dict{Participant, Dict{Reservoir, Float64}}(j => WeeklyAverageReservoirLevels[currentweek] for j in J)


Qadjs = create_Qadjs(R, K , 2)
Obligations, price, Qnoms = create_Obligations(T, J, R, Strategy, price_data, inflow_data, cuts, cutsOther, WaterCuts, WaterCutsOther, Initial_Reservoir, Initial_Individual_Reservoir, iteration_count_short, mu_up, mu_down, stage_count_bidding, scenario_count_prices, scenario_count_inflows, currentweek)
println("Obligations: $(Obligations)")
println("Price: $(price)")
println("Nominations: $(Qnoms)")
DisruptedNominations = Dict{Participant, Any}()
QnomOs = Dict{Participant, Any}()
QadjsFinal = Dict{Participant, Any}()
z_ups = Dict{Participant, Any}()
z_downs = Dict{Participant, Any}()
Revenues = Dict{Participant, Any}()
for j in J
    Qnom, QnomO, QadjFinal, z_up, z_down, Revenue = ShortTermProblemParametrized(R, K, J, j, Qnoms, Obligations, cuts, WaterCuts, Initial_Reservoir, Initial_Individual_Reservoir, mu_up, mu_down, price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, stage_count_short, T, 2)
    DisruptedNominations[j] = Qnom
    QnomOs[j] = QnomO
    QadjsFinal[j] = QadjFinal
    z_ups[j] = z_up
    z_downs[j] = z_down
    Revenues[j] = Revenue
end
"""
function ResultsToDataframeScheduling(savepath, J, R, Obligations, price, Qnoms, DisruptedNominations, currentweek)

    To save the results for later analysis, organize them inside a DataFrame.
    This is also helpful to do some statisical analysis, with functions from DataFrames.jl

"""
function ResultsToDataframeScheduling(savepath, J, R, Qadjs, Obligations, price, Revenues, QnomOs, DisruptedNominations, currentweek; save = true)
    column_names = ["week", "Participant", "Obligation", "price", ["Qnom_" * r.dischargepoint for r in R]..., ["QnomO_" * r.dischargepoint for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., "Revenue", "Upbalancing", "Downbalancing" ] 
    column_types = [Int64, String, Vector{Float64}, Vector{Float64}, [Float64 for r in R]..., [Float64 for r in R]...,  [Float64 for r in R]..., Float64, Float64, Float64]
    if isfile(savepath * "\\ShortTermScheduling\\DisruptedNominations.csv")
        println("DataFrame already exists. Add new Data....")
        df = CSV.File(savepath * "\\ShortTermScheduling\\DisruptedNominations.csv", types = column_types) |> DataFrame
        println(df)
    else
        println("Dataframe does not exist yet. Create DataFrame and fill with data...")
        df = DataFrame()
        for (name, type) in zip(column_names, column_types)
            df[!, name] = Vector{type}()
        end
        println(names(df))
    end
    for Qadj in Qadjs
        println(Qadj[R[1]])
        for j in J
            row = (
                week = currentweek,
                Participant = j.name,
                Obligation = Obligations[j],
                price = price,
                Qnom_Flasjon = DisruptedNominations[j][Qadj][(participant = j, reservoir = R[1])],
                Qnom_Holmsjon = DisruptedNominations[j][Qadj][(participant = j, reservoir = R[2])],
                QnomO_Flasjon = QnomOs[j][Qadj][j][R[1]],
                QnomO_Holmsjon = QnomOs[j][Qadj][j][R[2]],
                Qadj_Flasjon = Qadj[R[1]],
                Qadj_Holmsjon =  Qadj[R[2]],
                Revenue = Revenues[j][Qadj],
                Upbalancing = z_ups[j][Qadj],
                Downbalancing = z_downs[j][Qadj]               
            ) 
            push!(df, row)
        end
    end
    if save == true
        println("Results will be saved at $(savepath)...")
        CSV.write(savepath * "\\ShortTermScheduling\\DisruptedNominations.csv", df)
    end
    return df
end

df = ResultsToDataframeScheduling(savepath_results, J, R, Qadjs, Obligations, price, Revenues, QnomOs, DisruptedNominations, currentweek)