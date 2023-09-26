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

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)
R, K, J = read_data(filepath_Ljungan)
l_traj, f = AverageReservoirLevel(R, inflow_data)

const ColumnReservoir = Dict{Reservoir, String}(R[1] => "Flasjon Inflow", R[2] => "Holmsjon Inflow")
const stage_count_short = 2
const week = 2
const scenario_count_prices_medium = 3
const stage_count_medium = 52
const iteration_count_medium = 1500
const scenario_count_inflows_weekly = 1

NameToParticipant = Dict{String, Participant}(j.name => j for j in J)
PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)
Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j_loaded, MediumModelDictionary_O_loaded  = ReadMediumModel(savepath_watervalue, J, R, Ω_medium, P_medium, stage_count_medium, iteration_count_medium);

cuts_j = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, 2, stage_count_short) for j in J)
Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, week, stage_count_short) for j in J)
WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModelDictionary_j_loaded[j], cuts_j[j], week) for j in J)
WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModelDictionary_O_loaded[j], cutsOther[j], week) for j in J)


function SurfaceWaterValue(all_res::Vector{Reservoir}, V::SDDP.ValueFunction, j::Participant)
    R = collect(filter(r -> j.participationrate[r] > 0, all_res))
    if length(R) > 1
        bounds_1 = LinRange(0.0, R[1].currentvolume, 100)
        bounds_2 = LinRange(0.0, R[2].currentvolume, 100)
        f(x, y) = SDDP.evaluate(V, (Symbol("l[Holmsjon]") => y, Symbol("l[Flasjon]") => x))[1]
        println(f)
        surface(bounds_1, bounds_2, xlabel="$(R[1].dischargepoint)", ylabel="$(R[2].dischargepoint)", zlabel="Objective Value", f, title = "Water Value Function for $(j.name)", legend = false)
    else
        bounds_1 = LinRange(0.0, R[1].currentvolume, 100)
        g(x) = SDDP.evaluate(V, Dict(Symbol("l[Flasjon]") => x))[1]
        values_1 = g.(bounds_1)
        Plots.plot(bounds_1,values_1, legend = false, xlabel = "$(R[1].dischargepoint)", ylabel = "Objective Value", title = "Water Value Function for $(j.name)")
    end
end

"""
    PlotWaterValueCuts(all_res, V, j, cuts)

    Plot the linear cuts for the Water Value Functions. Cuts specifies the amount 
"""
function PlotWaterValueCuts(all_res::Vector{Reservoir}, V::SDDP.ValueFunction, j::Participant, cuts::Int64)
    R = collect(filter(r -> j.participationrate[r] > 0, all_res))
    ReservoirValues = Dict(r => collect(range(0, r.maxvolume, length=cuts)) for r in R)
    objvalues = [SDDP.evaluate(V, Dict(Symbol("l[$(r)]") => ReservoirValues[r][c] for r in R))[1] for c in 1:cuts]
    gradients = [SDDP.evaluate(V, Dict(Symbol("l[$(r)]") => ReservoirValues[r][c] for r in R))[2] for c in 1:cuts]
    Plots.plot(LinRange(0.0, R[1].maxvolume, 100), 
    [x -> objvalues[c] - min(objvalues...) -  sum(gradients[c][Symbol("l[$(r)]")] *(ReservoirValues[r][c] - x) for r in R) for c in 1:cuts],
    legend=false, 
    ylabel = "Objective Value",
    xlabel = L"Reservoir Volume $\left[\frac{m^3}{s}\right]$",
    title = "Water Value Cuts - $(j.name)",
    size = (1200, 800))
end
    
V_j  = Dict(j => SDDP.ValueFunction(MediumModelDictionary_j_loaded[j]; node = week) for j in J)
plot_WVC = PlotWaterValueCuts(R, V_j[j], j, 10)

j = J[3]
plot_WV = SurfaceWaterValue(R, V_j[j], j)
