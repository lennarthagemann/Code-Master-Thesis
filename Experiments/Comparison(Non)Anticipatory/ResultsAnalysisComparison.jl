
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