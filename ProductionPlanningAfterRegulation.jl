using JSON
using JuMP
using CPLEX
using Distributions
using LinearAlgebra
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation

# The nomination process is over, the flow has been adjusted by VF.
# What remains is to optimally allocate the water over data with minimal
# start-up costs and maximum revenue over the course of the day,
# While fulfilling adjusted flow constraints and fulfilling demand 
filepath_systemA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/Ljungan.json"

res, plants, parts = read_data(filepath_systemA)

j = parts[1]

a_1 = 0.5
a_2 = 1.5
AUCTION_LENGTH = 24
RES_UPPER_BOUND = 1000000
Qadj = 80
Qspill =  Dict(k => 105 for k in K)
P_swap = fill(18, AUCTION_LENGTH)
Pmin = Dict(k => 20 for k in K)
Pmax = Dict(k =>50 for k in K)
D = rand(Uniform(10, 100), AUCTION_LENGTH) 
mu_plus = rand(Uniform(5, 20), AUCTION_LENGTH)
mu_minus = rand(Uniform(2, 10), AUCTION_LENGTH)
Qinflow = 0
S = Dict(k => 1500 for k in K)
Res_ref = res[1].currentvolume

model = Model(CPLEX.Optimizer)
@variable(model, Qreal[1:AUCTION_LENGTH, K])
@variable(model, Preal[1:AUCTION_LENGTH, K])
@variable(model, z_plus[1:AUCTION_LENGTH] >= 0)
@variable(model, z_minus[1:AUCTION_LENGTH] >= 0)
@variable(model, u[1:AUCTION_LENGTH, K], Bin)
@variable(model, y[1:AUCTION_LENGTH, K], Bin) #linearize start-up costs@

@constraint(model, [k = K], Qadj == sum(Qreal[i,k] for i in 1:AUCTION_LENGTH)/AUCTION_LENGTH)
@constraint(model, [i = 1:AUCTION_LENGTH, k = K], a_1 * Qadj <= Qreal[i,k] <= a_2 * Qadj)
@constraint(model, [i = 1:AUCTION_LENGTH, k = K], u[i,k] * Pmin[k] <= Preal[i,k])
@constraint(model, [i = 1:AUCTION_LENGTH, k = K], Preal[i,k] <= u[i,k] * Pmax[k])
@constraint(model, [i = 1:AUCTION_LENGTH, k = K], Preal[i,k] <= Qreal[i,k] * k.equivalent)
@constraint(model, [i = 1:AUCTION_LENGTH, k = K], Preal[i,k] <= Qspill[k] * k.equivalent)
@constraint(model, [i = 1:AUCTION_LENGTH], D[i] == sum(Preal[i,k] for k in K) + P_swap[i] + (z_plus[i] - z_minus[i]))
for i in 1:AUCTION_LENGTH
    if i == 1
        @constraint(model, [k = K], y[i,k] >= u[i,k])
        @constraint(model, [k = K], y[i,k] >= 0)
    else
        @constraint(model, [k = K], y[i,k] >= u[i,k] - u[i-1,k])
        @constraint(model, [k = K], y[i,k] >= 0)
    end
end

@objective(model, Min, sum(sum(S[k] * y[:,k] for k in K)) + mu_plus'*z_plus + mu_minus'* z_minus)

optimize!(model)