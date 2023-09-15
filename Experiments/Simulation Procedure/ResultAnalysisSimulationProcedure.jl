#=
---------------------------------------------------------------------------------------------
In this file we analyze the Results from the Simulation Procedure

-Load in the DataFrame containing Results
-Plotting Function vor varying nomination and own nomination
-Plotting of how revenues change depending on adjsuted flow
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