"""
Test File to create new features for the water regulation process or the optimizaiton models
"""

using JuMP
using CPLEX
using Distributions
using LinearAlgebra
using Statistics
using Dates
using DataFrames
using SDDP
using Plots
import CSV
using JSON
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet(pwd() * "\\Water_Regulation\\WaterRegulation.jl")
using .WaterRegulation

filepath_Ljungan = pwd() * "\\Water_Regulation\\TestDataWaterRegulation\\Ljungan.json"
filepath_prices = pwd() *  "\\Data\\Spot Prices\\prices_df.csv"
filepath_results = pwd() * "\\Results\\LambdaZero\\"
res, plants, parts = read_data(filepath_Ljungan)

"""
Calculate cuts for generation function. Given a number of cuts, efficiency e and spill limit, return coefficients of the cuts 
by linear interpolation of the points.
"""

k = plants[1]
e = k.equivalent
Qspill = k.spill_reference_level
n = 3
function Generation_Cuts(Qspill::Float64, e::Float64, n::Int64)
    x = range(0,Qspill,n)
    effs  = range(e, e/2, n)
    y = x .* effs
    println("The generation curve points to interpolate are as follows: ")
    for i in 1:n
        println("(x = $(x[i]), y= $(y[i]))")
    end
    return x, y
end

kx, keffs = Generation_Cuts(Qspill, e, 5)