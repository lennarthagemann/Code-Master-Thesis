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


"""
function create_Qadjs(R, K, count)

    Helper Function to decide what kind of adjusted flow are to be used in the ShortTerScheduling Analysis.
    count gives the amount of different nominations to run through the short term scheduling
"""
function create_Qadjs(R::Vector{Reservoir}, K::Vector{HydropowerPlant}, count::Int64)::Vector{Dict{Reservoir, Float64}}
    Qadjs = Dict{Reservoir, Float64}[]
    for r in R
        max_flow = max([k.spillreference for k in filter(p -> p.reservoir in find_ds_reservoir(r) ,K)]...)
    end
    Qadj = Dict{Reservoir, Float64}(r => 0.0 for r in R)
    vals = [c * 1/counts * max_flow for c in 1:counts]
    for r in R
        Q[]
    end
    push!(Qadjs, Qadj)
    return Qadjs
end

"""
function ShortTermProblemParametrizedNomination(R, j, Obligation)

    The point of this function is to execute the ShortTermProblem for different input parameters. Mainly:

        -Different nomination by other people while Day Ahead Obligations are the same
    
    This leads to different power swap situations -> Technically this should be the same, but having to deliver a power swap might in some cases lead 
    to higher/lower nominations. The constraint of having to deliver a power swap at certaint times together with obligations
    might make it impossible to fulfill an obligation at a certain time.

    The Obligations will be taken from real obligations, resulting from either low/middle/high Bidding Curves.

    It will also be interesting to see if the participants are interested in th ereservoirs differently, even if the water volumes
    reaching their power plants are the same.
"""
function ShortTermProblemParametrized(R::Vector{Reservoir}, j::Participant, Obligation::Vector{Float64})
    

    Qadjs = create_Qadjs(R)

    QnomO1 = OthersNomination(Qnoms, Qadj1, J, R)
    Qnoms2 = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    Ω_NA_local, P_NA_local, _, _ = create_Ω_Nonanticipatory(price_data, inflow_data, scenario_count_prices, scenario_count_inflows, currentweek, R, stage_count_short)
    Qnom = ShortTermScheduling(R, j, J, Qref, Obligations[j], price, QnomO1[j], Ω_NA_local, P_NA_local, cuts[j], WaterCuts[j], iteration_count_short, mu_up,mu_down, T, stage_count_short) 
    for r in R
        if j.participationrate[r] > 0
            Qnoms2[(participant = j, reservoir = r)] = Qnom[r]
        else
            Qnoms2[(participant = j, reservoir = r)] = Qref[r]
        end
    end
    
    
    return
end