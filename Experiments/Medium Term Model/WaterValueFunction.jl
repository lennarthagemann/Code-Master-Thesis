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
const stage_count_medium = 52
const stage_count_short = 2
const scenario_count_prices_medium = 3
const scenario_count_inflows_weekly = 3
const iterations_medium = 100
const week = 2

price_data = prepare_pricedata(filepath_prices)
inflow_data = prepare_inflowdata(filepath_inflows)
l_traj, f = AverageReservoirLevel(R, inflow_data)

PriceScenariosMedium = Price_Scenarios_Medium(price_data, scenario_count_prices_medium)
InflowScenariosMedium = Inflow_Scenarios_Medium(inflow_data, ColumnReservoir, scenario_count_inflows_weekly, R)

Ω_medium, P_medium =  create_Ω_medium(PriceScenariosMedium, InflowScenariosMedium, R);
MediumModelDictionary_j, MediumModelDictionary_O = MediumModelsAllParticipants(J, R, Ω_medium, P_medium, stage_count_medium, iterations_medium; print_level = 1)
for j in J
    SDDP.write_cuts_to_file(MediumModelDictionary_j[j], savepath_watervalue * "\\Participant\\$(j).json")
    SDDP.write_cuts_to_file(MediumModelDictionary_O[j], savepath_watervalue * "\\OtherParticipant\\$(j).json")
end
# SaveMediumModel(savepath_watervalue * "\\Participant.jls", MediumModelDictionary_j)
# SaveMediumModel(savepath_watervalue * "\\OtherParticipant.jls", MediumModelDictionary_O)

cuts_j = Dict(j => ReservoirLevelCuts(R, j.plants, j, f, week, 7) for j in J)
Others = Dict(j => OtherParticipant(J,j,R)[1] for j in J)
cutsOther = Dict(j => ReservoirLevelCuts(R, Others[j].plants, Others[j], f, week, stage_count_short) for j in J)
WaterCuts = Dict(j => WaterValueCuts(R, j, MediumModelDictionary_j[j], cuts_j[j], week) for j in J)
WaterCutsOther = Dict(j => WaterValueCuts(R, Others[j], MediumModelDictionary_O[j], cutsOther[j], week) for j in J)


V_j = Dict(j => SDDP.ValueFunction(MediumModelDictionary_j[j]; node = 52) for j in J)
bounds_1 = LinRange(0.0, R[1].currentvolume, 100)
bounds_2 = LinRange(0.0, R[2].currentvolume, 100)
WV(x, y) = SDDP.evaluate(V_j[J[2]], Dict(Symbol("l[Holsmjon]") => y, Symbol("l[Flasjon]") => x))[1]
surface(bounds_1, bounds_2, WV, legend = false)

g(x) = SDDP.evaluate(V_j, (Symbol("l[Holsmjon]") => x, Symbol("l[Flasjon]") => 0))[1]
h(x) = SDDP.evaluate(V_j, (Symbol("l[Holsmjon]") => 0, Symbol("l[Flasjon]") => x))[1]

Plots.plot([bounds_1, bounds_2], [g, h],
    xlabel = "Reservoir Level", ylabel = "Objective Value", legend=false, layout = (1, 2), size= (1000,600), show=true)

cuts = 5
ReservoirValues = Dict(r => collect(range(0, r.maxvolume, length=cuts))  for r in R)
objvalues = [SDDP.evaluate(V_j, Dict(Symbol("l[$(r)]") => ReservoirValues[r][c] for r in R))[1] for c in 1:cuts]
gradients = [SDDP.evaluate(V_j, Dict(Symbol("l[$(r)]") => ReservoirValues[r][c] for r in R))[2] for c in 1:cuts]

Plots.plot(LinRange(0.0, R[1].maxvolume, 100), 
    [x -> objvalues[c] - min(objvalues...) -  gradients[c][Symbol("l[$(R[1])]")] *(ReservoirValues[R[1]][c] - x) - gradients[c][Symbol("l[$(R[2])]")] * (ReservoirValues[R[2]][c] - x) for c in 1:cuts],
    legend=false, 
    ylabel = "Objective Value",
    xlabel = "Reservoir Level",
    title = "Water Value Cuts")

Plots.plot(LinRange(0.0, R[1].maxvolume, 100), 
    [x ->  WaterCuts[j][c].e1 - min([WaterCuts[j][cut].e1 for cut in cuts_j[j]]...) - sum(WaterCuts[j][c].e2[Symbol("l[$(r)]")] *(c[r] - x) for r in R) for c in cuts_j[j]],
    legend=false, 
    ylabel = "Objective Value",
    xlabel = "Reservoir Level",
    title = "Water Value Cuts")


"""
    SurfaceWaterValue(all_res, V, j)

    Plot the Objective Value of the Medium Term Model for the two reservoirs as surface plot.
"""
function SurfaceWaterValue(all_res::Vector{Reservoir}, V::SDDP.ValueFunction, j::Participant)
    R = collect(filter(r -> j.participationrate[r] > 0, all_res))
    if length(R) > 1
        bounds_1 = LinRange(0.0, R[1].currentvolume, 100)
        bounds_2 = LinRange(0.0, R[2].currentvolume, 100)
        f(x, y) = SDDP.evaluate(V, (Symbol("l[Holsmjon]") => y, Symbol("l[Flasjon]") => x))[1]
        surface(bounds_1, bounds_2, f, legend = false)
    else
        bounds_1 = LinRange(0.0, R[1].curremtvolume, 100)
        f(x) = SDDP.evaluate(V, (Symbol("l[Flasjon]") => x))[1]
        surface(bounds_1, f, legend = false)
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
    xlabel = "Reservoir Level",
    title = "Water Value Cuts")
end

