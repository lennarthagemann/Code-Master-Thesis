
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
using Statistics
using Dates
using DataFrames
using CSV
using JSON
using PlotlyJS
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end

includet(pwd() * "\\Water_Regulation\\WaterRegulation.jl")
using .WaterRegulation

# Custom function to parse the string into a Vector{Float64}
function parse_vector_string(s::AbstractString)
    values_str = replace(s, r"[\[\]]" => "")
    values = split(values_str, ",")

    # Convert the string elements to Float64 and create a Vector{Float64}
    vector_obj = parse.(Float64, values)
    return vector_obj
end

Base.parse(::Type{Vector{Float64}}, s::AbstractString) = parse_vector_string(s)
Base.tryparse(::Type{Vector{Float64}}, s::AbstractString) = parse_vector_string(s)


const filepath_Ljungan = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
const savepath_experiment_nominations = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\NonanticipatoryVsAnticipatory\\NominationsBounded.csv"
const savepath_experiment_obligations = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\NonanticipatoryVsAnticipatory\\ObligationsBounded.csv"
const savepath_experiment_reservoirs = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\NonanticipatoryVsAnticipatory\\ReservoirsBounded.csv"


R, K, J = read_data(filepath_Ljungan)

function LoadResultsDataFrameNominations(R::Vector{Reservoir}, J::Vector{Participant}, savepath_experiment::String)
    column_names_df_nominations = ["week", ["Strategy_" * j.name for j in J]..., ["Qnom1_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qnom2_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., ["P_Swap_" * j.name * "_" * r.dischargepoint for j in J for r in R]...]
    column_types_df_nominations = [Int64, [String for j in J]..., [Float64 for j in J for r in R]..., [Float64 for j in J for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]...]
   df = CSV.File(savepath_experiment, types = column_types_df_nominations) |> DataFrame
    return df
end

function LoadResultsDataFrameObligations(R::Vector{Reservoir}, J::Vector{Participant}, savepath_experiment::String)
    column_names_df_Obligations = ["week", ["Strategy_" * j.name for j in J]..., ["Obligations_" * j.name for j in J]..., ["z_up_" * j.name for j in J]..., ["z_down_" *j.name for j in J]..., ["Revenue_"*j.name for j in J]...]
    column_types_df_Obligations = [Int64, [String for j in J]..., [Vector{Float64} for j in J]..., [Vector{Float64} for j in J]..., [Vector{Float64} for j in J]..., [Float64 for j in J]...]
    df = CSV.File(savepath_experiment, types = column_types_df_Obligations) |> DataFrame
    return df
end
function LoadResultsDataFrameReservoirs(R::Vector{Reservoir}, J::Vector{Participant}, savepath_experiment::String)
    column_names_df_Reservoirs = ["week", ["Strategy_" * j.name for j in J]...,  ["l_real_" * r.dischargepoint for r in R]..., ["l_ind_" * j.name * "_" * r.dischargepoint for j in J for r in R]...]
    column_types_df_Reservoirs = [Int64, [String for j in J]...,  [Float64 for r in R]..., [Float64 for j in J for r in R]...,]
       df = CSV.File(savepath_experiment, types = column_types_df_Reservoirs) |> DataFrame
    return df
end

results_df_nominations = LoadResultsDataFrameNominations(R, J, savepath_experiment_nominations)
results_df_obligations = LoadResultsDataFrameObligations(R, J, savepath_experiment_obligations)
results_df_reservoirs = LoadResultsDataFrameReservoirs(R, J, savepath_experiment_reservoirs)

function combineDF(J::Vector{Participant}, R::Vector{Reservoir}, results_df_nominations, results_df_obligations, results_df_reservoirs)
    results_df = select(hcat(results_df_nominations, results_df_obligations,  results_df_reservoirs, makeunique=true), [:week, :Strategy_Sydkraft, :Strategy_Fortum, :Strategy_Statkraft,
    :Qnom2_Sydkraft_Flasjon, :Qnom2_Fortum_Flasjon, :Qnom2_Fortum_Holmsjon, :Qnom2_Statkraft_Flasjon, :Qnom2_Statkraft_Holmsjon, :Revenue_Sydkraft, :Revenue_Fortum, :Revenue_Statkraft,
    ])
    results_df = select(results_df, :week, :Strategy_Sydkraft, :Strategy_Fortum, :Strategy_Statkraft,
    :Qnom2_Sydkraft_Flasjon, :Qnom2_Fortum_Flasjon, :Qnom2_Fortum_Holmsjon, :Qnom2_Statkraft_Flasjon, :Qnom2_Statkraft_Holmsjon, :Revenue_Sydkraft, :Revenue_Fortum, :Revenue_Statkraft,
    [:Qnom2_Sydkraft_Flasjon, :Revenue_Sydkraft,] => ((x1, x2) -> x2 ./ (x1 * J[1].participationrate[R[1]])) => :RevenueScaled_Sydkraft,
    [:Qnom2_Fortum_Flasjon, :Qnom2_Fortum_Holmsjon, :Revenue_Sydkraft,] => ((x1, x2, x3) -> x3 ./ ((x1 .+ x2) *J[2].participationrate[R[1]])) => :RevenueScaled_Fortum,
    [:Qnom2_Statkraft_Flasjon, :Qnom2_Statkraft_Holmsjon, :Revenue_Sydkraft,] => ((x1, x2, x3) -> x3 ./ ((x1 .+ x2) * J[3].participationrate[R[1]])) => :RevenueScaled_Statkraft, 
    )
    ScaledRevenuesNonanticipatory = Dict{Participant, DataFrame}() 
    ScaledRevenuesAnticipatory = Dict{Participant, DataFrame}() 
    ScaledRevenuesNonanticipatory[J[1]]  = select(subset(results_df, :Strategy_Sydkraft  => x -> x.== "Nonanticipatory"), :RevenueScaled_Sydkraft)
    ScaledRevenuesAnticipatory[J[1]] = select(subset(results_df, :Strategy_Sydkraft  => x -> x.== "Anticipatory", :RevenueScaled_Sydkraft => x -> x .> -1e6), :RevenueScaled_Sydkraft)
    ScaledRevenuesNonanticipatory[J[2]] = select(subset(results_df, :Strategy_Fortum  => x -> x.== "Nonanticipatory"), :RevenueScaled_Fortum)
    ScaledRevenuesAnticipatory[J[2]] = select(subset(results_df, :Strategy_Fortum  => x -> x.== "Anticipatory"), :RevenueScaled_Fortum)
    ScaledRevenuesNonanticipatory[J[3]] = select(subset(results_df, :Strategy_Statkraft  => x -> x.== "Nonanticipatory"), :RevenueScaled_Statkraft)
    ScaledRevenuesAnticipatory[J[3]] = select(subset(results_df, :Strategy_Statkraft  => x -> x.== "Anticipatory"), :RevenueScaled_Statkraft)
    return results_df, ScaledRevenuesNonanticipatory, ScaledRevenuesAnticipatory 
end

results_df, ScaledRevenuesNonanticipatory, ScaledRevenuesAnticipatory = combineDF(J, R, results_df_nominations, results_df_obligations, results_df_reservoirs)

plot([
    histogram(ScaledRevenuesAnticipatory[J[1]], x=:RevenueScaled_Sydkraft, opacity=0.9, nbinsx = 10, name = J[1].name),
    histogram(ScaledRevenuesAnticipatory[J[2]], x=:RevenueScaled_Fortum, opacity=0.9, nbinsx = 10, name = J[2].name),
    histogram(ScaledRevenuesAnticipatory[J[3]], x=:RevenueScaled_Statkraft, opacity=0.9, nbinsx = 10, name = J[3].name)],
    Layout(title = "Statistical Distribution of Revenues - Anticipatory", xaxis_title= "Revenue", xaxis_range = [-2000, 4500], yaxis_title = "Count" ))

plot([
    histogram(ScaledRevenuesNonanticipatory[J[1]], x=:RevenueScaled_Sydkraft, opacity=0.9, nbinsx = 10, name = J[1].name),
    histogram(ScaledRevenuesNonanticipatory[J[2]], x=:RevenueScaled_Fortum, opacity=0.9, nbinsx = 10, name = J[2].name),
    histogram(ScaledRevenuesNonanticipatory[J[3]], x=:RevenueScaled_Statkraft, opacity=0.9, nbinsx = 10, name = J[3].name)],
    Layout(title = "Statistical Distribution of Revenues - Nonanticipatory", xaxis_title= "Revenue", xaxis_range = [-2000, 4500], yaxis_title = "Count" ))