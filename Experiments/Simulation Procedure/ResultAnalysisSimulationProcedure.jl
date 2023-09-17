#=
---------------------------------------------------------------------------------------------
In this file we analyze the Results from the Simulation Procedure

-Load in the DataFrame containing Results
-Plotting Function vor varying nomination and own nomination
-Plotting of how revenues change depending on adjsuted flow
---------------------------------------------------------------------------------------------
=#


using Distributions
using Statistics
using Dates
using DataFrames
using CSV
using JSON
using PlotlyJS
using LaTeXStrings
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
const savepath_experiment = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\SingleVsIndividual\\SingleVsIndividualBounded.csv"

R, K, J = read_data(filepath_Ljungan)

function LoadResultsDataFrame(R::Vector{Reservoir}, J::Vector{Participant}, savepath_experiment::String)
    column_names = ["week", "day", ["Strategy_" * j.name for j in J]..., ["Reservoir_Single_" * r.dischargepoint for r in R]..., ["Reservoir_VF_" * r.dischargepoint for r in R]..., ["Individual_Reservoir_" * j.name * "_" * r.dischargepoint for j in J for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., ["Q_single_" * r.dischargepoint for r in R]..., ["Qnom_$(j.name)_$(r.dischargepoint)" for j in J for r in R]..., ["P_Swaps_$(j.name)_$(r.dischargepoint)" for j in J for r in R]..., ["Obligations_" * j.name for j in J]...,"Obligations_single", ["VF_Revenues_" * j.name for j in J]..., "Single_Revenue", ["Split_Revenues_" * j.name for j in J]..., ["z_ups_" * j.name for j in J]..., ["z_downs_" * j.name for j in J]... ]
    column_types = [Int64, Int64, [String for j in J]..., [Float64 for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]..., [Float64 for r in R]..., [Float64 for r in R]..., [Float64 for j in J for r in R]..., [Float64 for j in J for r in R]..., [Vector{Float64} for j in J]..., Vector{Float64},  [Float64 for j in J]..., Float64, [Float64 for j in J]..., [Float64 for j in J]..., [Float64 for j in J]...]
    df = CSV.File(savepath_experiment, types = column_types) |> DataFrame
    return df
end

function ScaledRevenues(results_df::DataFrame, K::Vector{HydropowerPlant}, J::Vector{Participant}, R::Vector{Reservoir})
    revenues_df_Nonanticipatory = select(subset(results_df, :Strategy_Sydkraft => x -> x .=="Nonanticipatory"), :Q_single_Flasjon, :Q_single_Holmsjon, :Qnom_Sydkraft_Flasjon,  :Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon, :Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon,  :Single_Revenue, :VF_Revenues_Sydkraft, :VF_Revenues_Fortum, :VF_Revenues_Statkraft, :Split_Revenues_Sydkraft, :Split_Revenues_Fortum, :Split_Revenues_Statkraft)
    revenues_df_Anticipatory = select(subset(results_df, :Strategy_Sydkraft => x -> x .=="Anticipatory"), :Q_single_Flasjon, :Q_single_Holmsjon, :Qnom_Sydkraft_Flasjon,  :Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon, :Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon, :Single_Revenue, :VF_Revenues_Sydkraft, :VF_Revenues_Fortum, :VF_Revenues_Statkraft, :Split_Revenues_Sydkraft, :Split_Revenues_Fortum, :Split_Revenues_Statkraft)

    revenues_df_Nonanticipatory = select(revenues_df_Nonanticipatory, :Single_Revenue, :VF_Revenues_Sydkraft, :VF_Revenues_Fortum, :VF_Revenues_Statkraft, :Split_Revenues_Sydkraft, :Split_Revenues_Fortum, :Split_Revenues_Statkraft, :Q_single_Flasjon, :Q_single_Holmsjon, :Qnom_Sydkraft_Flasjon,  :Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon, :Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon,
    [:Single_Revenue, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ (x2 .* J[1].participationrate[R[1]] + x3 .* (J[2].participationrate[R[1]] + J[3].participationrate[R[1]]))) =>  :Scaled_Revenue,
    [:VF_Revenues_Sydkraft ,:Qnom_Sydkraft_Flasjon] => ((x1,x2) -> min.(x1 ./ (x2 * J[1].participationrate[R[1]]), 4000)) =>  :Scaled_Revenue_Sydkraft,
    [:VF_Revenues_Fortum,:Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[2].participationrate[R[1]])) =>  :Scaled_Revenue_Fortum,
    [:VF_Revenues_Statkraft,:Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[3].participationrate[R[1]])) =>  :Scaled_Revenue_Statkraft,

    [:Split_Revenues_Sydkraft ,:Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[1].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Sydkraft,
    [:Split_Revenues_Fortum, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[2].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Fortum,
    [:Split_Revenues_Statkraft, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[3].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Statkraft,                                                                               
    )    
    println(select( revenues_df_Nonanticipatory))
    revenues_df_Anticipatory = select(revenues_df_Anticipatory, :Single_Revenue, :VF_Revenues_Sydkraft, :VF_Revenues_Fortum, :VF_Revenues_Statkraft, :Split_Revenues_Sydkraft, :Split_Revenues_Fortum, :Split_Revenues_Statkraft, :Q_single_Flasjon, :Q_single_Holmsjon, :Qnom_Sydkraft_Flasjon,  :Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon, :Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon,
    [:Single_Revenue, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ (x2 .* J[1].participationrate[R[1]] + x3 .* (J[2].participationrate[R[1]] + J[3].participationrate[R[1]]))) =>  :Scaled_Revenue,
    
    [:VF_Revenues_Sydkraft ,:Qnom_Sydkraft_Flasjon] => ((x1,x2) -> x1 ./ (x2 * J[1].participationrate[R[1]])) =>  :Scaled_Revenue_Sydkraft,
    [:VF_Revenues_Fortum, :Qnom_Fortum_Flasjon, :Qnom_Fortum_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[2].participationrate[R[1]])) =>  :Scaled_Revenue_Fortum,
    [:VF_Revenues_Statkraft, :Qnom_Statkraft_Flasjon, :Qnom_Statkraft_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[3].participationrate[R[1]])) =>  :Scaled_Revenue_Statkraft,
    
    [:Split_Revenues_Sydkraft, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[1].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Sydkraft,
    [:Split_Revenues_Fortum, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[2].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Fortum,
    [:Split_Revenues_Statkraft, :Q_single_Flasjon, :Q_single_Holmsjon] => ((x1,x2,x3) -> x1 ./ ((x2 .+ x3) * J[3].participationrate[R[1]])) =>  :Scaled_Split_Revenue_Statkraft,                                                                               
    )    
    
    revenues_df_Nonanticipatory = select(subset(revenues_df_Nonanticipatory, :Scaled_Revenue_Sydkraft =>  x -> x .> -1e6),  [:Scaled_Revenue, :Scaled_Revenue_Sydkraft, :Scaled_Revenue_Fortum, :Scaled_Revenue_Statkraft, :Scaled_Split_Revenue_Sydkraft, :Scaled_Split_Revenue_Fortum, :Scaled_Split_Revenue_Statkraft])
    revenues_df_Anticipatory = select(subset(revenues_df_Anticipatory, :Scaled_Revenue_Sydkraft => x  -> x .> -1e6 ), [:Scaled_Revenue, :Scaled_Revenue_Sydkraft, :Scaled_Revenue_Fortum, :Scaled_Revenue_Statkraft,  :Scaled_Split_Revenue_Sydkraft, :Scaled_Split_Revenue_Fortum, :Scaled_Split_Revenue_Statkraft])
    return revenues_df_Nonanticipatory, revenues_df_Anticipatory
end

results_df = LoadResultsDataFrame(R, J, savepath_experiment)
revenues_df_Nonanticipatory, revenues_df_Anticipatory  = ScaledRevenues(results_df, K, J, R)

plot([
    histogram(revenues_df_Nonanticipatory, x=:Scaled_Revenue, opacity=0.9, name = "Single Participant"),
    histogram(revenues_df_Nonanticipatory, x=:Scaled_Revenue_Sydkraft, opacity=0.9, name = J[1].name),
    histogram(revenues_df_Nonanticipatory, x=:Scaled_Revenue_Fortum, opacity=0.9, name = J[2].name),
    histogram(revenues_df_Nonanticipatory, x=:Scaled_Revenue_Statkraft, opacity=0.9, name = J[3].name)],
    Layout(title = L"Statistical Distribution of Revenues per \$\\frac{m^3}{s}\$ - Nonanticipatory", xaxis_title= "Revenue", yaxis_title = "Count" ))

plot([
    histogram(revenues_df_Anticipatory, x=:Scaled_Revenue, opacity=0.9, name = "Single Participant"),
    histogram(revenues_df_Anticipatory, x=:Scaled_Revenue_Sydkraft, opacity=0.9, name = J[1].name),
    histogram(revenues_df_Anticipatory, x=:Scaled_Revenue_Fortum, opacity=0.9, name = J[2].name),
    histogram(revenues_df_Anticipatory, x=:Scaled_Revenue_Statkraft, opacity=0.9, name = J[3].name)],
    Layout(title = "Statistical Distribution of Revenues per \\frac{m^3}{s}  - Anticipatory", xaxis_title= "Revenue", yaxis_title = "Count" ))


select(results_df, :Single_Revenue , [:Q_single_Flasjon, :Q_single_Holmsjon] => ((x1, x2) ->  x1 .+ x2) => :Qsingle, [:Qadj_Flasjon, :Qadj_Holmsjon] => ((x1, x2) -> x1 .+ x2) => :Qadj, [:VF_Revenues_Sydkraft, :VF_Revenues_Fortum, :VF_Revenues_Statkraft] => ((x1, x2, x3) ->  x1 .+ x2 .+x3) => :VF_Revenues)
select(revenues_df_Nonanticipatory, r"Scaled")