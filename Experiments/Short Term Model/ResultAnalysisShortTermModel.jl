
#=
---------------------------------------------------------------------------------------------
In this file we analyze the Results from the ShortTermModel Analysis

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
const savepath_experiment = "C:\\Users\\lenna\\OneDrive - NTNU\\Code Master Thesis\\Experiments\\Results\\ShortTermScheduling\\DisruptedNominationsSimple.csv"

function LoadResultsDataFrame(R::Vector{Reservoir}, savepath_experiment)
    column_names = ["week", "Participant", "Obligation", "price", ["Qnom_" * r.dischargepoint for r in R]..., ["QnomO_" * r.dischargepoint for r in R]..., ["Qadj_" * r.dischargepoint for r in R]..., "Revenue", "Upbalancing", "Downbalancing" ] 
    column_types = [Int64, String, Vector{Float64}, Vector{Float64}, [Float64 for r in R]..., [Float64 for r in R]...,  [Float64 for r in R]..., Float64, Float64, Float64]
    df = CSV.File(savepath_experiment, types = column_types) |> DataFrame
    return df
end

R, K, J = read_data(filepath_Ljungan)
results_df = LoadResultsDataFrame(R, savepath_experiment)

j = J[2]
week = 2

subdf = filter(row -> (row.week == 2 && row.Participant == j.name), select(results_df, ["Participant", "week", "Obligation", "Qnom_Flasjon", "Qnom_Holmsjon", "Qadj_Flasjon", "Qadj_Holmsjon", "QnomO_Flasjon", "QnomO_Holmsjon"]))
function PlotScatterOneReservoir(j::Participant)
    subdf = filter(row -> (row.week == 2 && row.Participant == j.name), select(results_df, ["Participant", "week", "Obligation", "Qnom_Flasjon", "Qnom_Holmsjon", "Qadj_Flasjon", "Qadj_Holmsjon", "QnomO_Flasjon", "QnomO_Holmsjon"]))
    p = plot([
        bar(subdf, x= nrow(subdf), y=:Qnom_Flasjon, kind="scatter", mode="markers", marker_color="blue", name="Own Nomination" ),
        bar(subdf, x= nrow(subdf), y=:QnomO_Flasjon, kind="scatter", mode="markers", marker_color="grey", name="Other's Nomination" ),
        bar(subdf, x= nrow(subdf), y=:Qadj_Flasjon, kind="scatter", mode="markers", marker_color="lightblue", name="Adjusted Flow" ),
    ], Layout(xaxis_tile = "", yaxis_title = "Discharge Flasjon", title="Nominations from Short Term Scheduling - Sydkraft"))
    return p
end
# p = PlotScatterOneReservoir(J[1])

function PlotScatterTwoReservoir(j::Participant, week::Int64)
    subdf = filter(row -> (row.week == week && row.Participant == j.name), select(results_df, ["Participant", "week", "Obligation", "Qnom_Flasjon", "Qnom_Holmsjon", "Qadj_Flasjon", "Qadj_Holmsjon", "QnomO_Flasjon", "QnomO_Holmsjon"]))
    trace1 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:Qnom_Flasjon, z=:Qnom_Holmsjon, kind="scatter3d", mode="markers", marker_color="blue", name="Own Nomination")
    trace2 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:QnomO_Flasjon, z=:QnomO_Holmsjon, kind="scatter3d", mode="markers", marker_color="grey", name="Other's Nomination")
    trace3 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:Qadj_Flasjon, z=:Qadj_Holmsjon, kind="scatter3d", mode="markers", marker_color="lightblue", name="Adjusted Flow")
    Traces = [trace1, trace2, trace3]
    p3d = plot(Traces, Layout(scene = attr(
        xaxis_title="",
        yaxis_title="Nomination Flasjon",
        zaxis_title="Nomination Holmsjon"), title="Nominations from Short Term Scheduling"))
    return p3d
end
# p3d = PlotScatterTwoReservoir(J[2], 20)

function PlotScatterTwoReservoir(j::Participant, week::Int64)
    subdf = filter(row -> (row.week == week && row.Participant == j.name), select(results_df, ["Participant", "week", "Obligation", "Qnom_Flasjon", "Qnom_Holmsjon", "Qadj_Flasjon", "Qadj_Holmsjon", "QnomO_Flasjon", "QnomO_Holmsjon"]))
    trace1 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:Qnom_Flasjon, z=:Qnom_Holmsjon, kind="scatter3d", mode="markers", marker_color="blue", name="Own Nomination")
    trace2 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:QnomO_Flasjon, z=:QnomO_Holmsjon, kind="scatter3d", mode="markers", marker_color="grey", name="Other's Nomination")
    trace3 = PlotlyJS.scatter(subdf, x= 1:nrow(subdf), y=:Qadj_Flasjon, z=:Qadj_Holmsjon, kind="scatter3d", mode="markers", marker_color="lightblue", name="Adjusted Flow")
    Traces = [trace1, trace2, trace3]
    p3d = plot(Traces, Layout(scene = attr(
        xaxis_title="",
        yaxis_title="Nomination Flasjon",
        zaxis_title="Nomination Holmsjon"), title="Nominations from Short Term Scheduling"))
    return p3d
end
function Revenue3dScatter(j::Participant, week::Int64, results_df::DataFrame)
    subdf = filter(row -> (row.week == week && row.Participant == j.name), select(results_df, :Participant, :week, :Revenue, :Qnom_Flasjon, :Qnom_Holmsjon, :Qadj_Flasjon, :Qadj_Holmsjon, :QnomO_Flasjon, :QnomO_Holmsjon))
    transform!(subdf, [:Qnom_Flasjon, :Qnom_Holmsjon] => (+) => :Qnom)
    transform!(subdf, [:Qadj_Flasjon, :Qadj_Holmsjon] => (+) => :Qadj)
    transform!(subdf, [:QnomO_Flasjon, :QnomO_Holmsjon] => (+) => :QnomO)
    select!(subdf, :Participant, :week, :Revenue, :Qnom, :Qadj, :QnomO)
    trace1 = PlotlyJS.scatter(subdf, x=:Revenue, y=:Qnom, z=:QnomO, kind="scatter3d", mode="markers", marker_color="blue", name="Own Nomination")
    Traces = [trace1]
    p3d = plot(Traces, Layout(scene = attr(
        xaxis_title="Revenue",
        yaxis_title="Nomination Flasjon",
        zaxis_title="Nomination Holmsjon"), title="Revenue under Disruptions"))
    return p3d
end 

Revenue3dScatter(j, 20, results_df)