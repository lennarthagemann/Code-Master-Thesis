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

"""
_____________________________________________________________________________________
-                                   Parameters                                      -
_____________________________________________________________________________________

"""

R, K, J = read_data(filepath_Ljungan)
println("Reservoirs: ", R)
println("Power Plants: ", K)
println("Participants: ",J) 
stages = 2
T = 24
PPoints = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
I = length(PPoints) - 1
Qref = Dict{Reservoir, Float64}(r => 10.0 for r in R)
scenario_count = 10
Prices = [floor.(rand(T), sigdigits=3) for i in 1:scenario_count]
Inflows = [10.0]
# StartUp Costs
S = 0.1
# Cost for Up and Down Balancing
mu_up = 0.7
mu_down = 0.3


#=
_____________________________________________________________________________________
-                                Helper Functions                                   -
_____________________________________________________________________________________

=#


#=
Create Uncertainty Set(s) for the nonanticipatory problem respectively.
For the scope of this thesis we constrain ourselves to the cartesian products of inflow and prices and don't include any dependencies of both values.
=#

"""
From Inflow and Price Scenarios, create the Uncertainty Set(s) which make up the scenario tree. 
"""
function Uncertainty_Nonanticipatory(Prices, Inflows)
    Omega = [(price = p, inflow = v) for p in Prices for v in Inflows]
    P = [1/length(Omega) for om in Omega]
    return Omega, P
end
"""
From Inflow and Price Scenarios, create the Uncertainty Set(s) which make up the scenario tree. 
For Inflow and Price Scenarios, think about how the other 
"""

function Uncertainty_Anticipatory(Prices, Inflows)
    Qnom_O = Dict{Reservoir, Float64}(r => 0.3 for r in R) # Here the nonanticipatory Bidding Problem has to be solved.
    Omega = [(price = p, inflow = v, nomination = Qnom_O) for p in Prices for v in Inflows]
    P = [1/length(Omega) for om in Omega]
    return Omega, P
end

"""
Simulate Market Clearing: For a Price Realization and given Bidding Curves, determine the obligation for each producer.
    1. Find the price segment where the realized price lies
    2. Select the bids from surrounding Price Points of segment
    3. Obtain delivery through linear interpolation
"""
function MarketClearing(price::Array{Float64}, PPoints::Array{Float64}, BidCurves::Dict{Participant, Dict{Int64, Vector{Float64}}}, J::Array{Participant}; T::Int64 = T)::Dict{Participant, Dict{Int64, Float64}}
    @assert(length(PPoints) > 1)
    CP = Dict{Int64, Int64}(t => 0 for t in 1:T)
    for t in 1:T
        for i in 1:length(PPoints)-1
            if price[t] >= PPoints[i] && price[t] <= PPoints[i+1]
                CP[t] = i
                break
            end
        end
    end

    y = Dict{Participant, Dict{Int64, Float64}}(j => Dict(t => 0.0 for t in 1:T) for j in J)
    for j in J
        for t in 1:T
            y[j][t] =(price[t] - PPoints[CP[t]])/(PPoints[CP[t]+1] - PPoints[CP[t]]) * BidCurves[j][t][CP[t]] + (PPoints[CP[t]+1] - price[t])/(PPoints[CP[t]+1] - PPoints[CP[t]]) * BidCurves[j][t][CP[t]+1]
        end
    end
    return y
end
"""
After A round of nominations, determine the aggregated others nomination from the adjusted flow and own nomination.

"""
function OtherNominations(J::Array{Participant}, R::Array{Reservoir}, Nominations::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qadj::Dict{Reservoir, Float64})
    QnomO = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}((participant = j, reservoir = r) => 0.0 for j in J for r in R)
    for j in J
        O, K_O = OtherParticipant(J, j ,R)
        for r in filter(r -> j.participationrate[r] > 0.0, R)
            QnomO[(participant = j, reservoir = r)] = (Qadj[r] * (j.participationrate[r] + O.participationrate[r]) - Nominations[(participant = j, reservoir = r)] * j.participationrate[r])/O.participationrate[r]
        end
    end
    return QnomO
end

"""
After all Optimziation Problems have been completed, it is time to calculate how much money each producer has made from that day.
For a fixed price and delivery, calculate how much money has been made through sales on the spot market. To this add the up and downregulation on the balancing market.
"""
function CalculateRevenues(
    J::Vector{Participant},
    y::Dict{Participant, Dict{Int64, Float64}},
    price::Vector{Float64},
    Upregulation::Dict{Participant, Vector{Float64}},
    Downregulation::Dict{Participant, Vector{Float64}};
    mu_up = 1.0,
    mu_down = 0.1
    )::Dict{Participant, Float64}
    Revenues = Dict{Participant, Float64}(j => 0.0 for j in J)
    for j in J
        Revenues[j] = sum(price[t] * y[j][t] + Downregulation[j][t] * mu_down - Upregulation[j][t] * mu_up for t in 1:T)
    end
    return Revenues
end

#=
_____________________________________________________________________________________
-                                   Simulations                                     -
_____________________________________________________________________________________

=#

Omega, P = Uncertainty_Anticipatory(Prices, Inflows)

# --------------------------- Bidding and Nominaionts -------------------------------

Nominations = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}((participant = j, reservoir = r) => 0.0 for r in R for j in J)
BidCurves = Dict{Participant, Dict{Int64, Vector{Float64}}}(j => Dict(t => [0.0] for t in 1:T) for j in J)
for j in J
    Qnom_temp, BidCurve_temp = Nonanticipatory_Bidding(R, j, PPoints, Omega, P)
    for r in filter(r -> j.participationrate[r] > 0, R)
        Nominations[(participant = j, reservoir = r)] = Qnom_temp[r]
    end
    BidCurves[j] = BidCurve_temp
end

# Plots.plot([[BidCurves[j][i] for j in J] for i in 1:T], layout = 3)

# Sample Price and market clearing
sample = Prices[rand(1:length(Prices))]
y = MarketClearing(sample, PPoints, BidCurves, J)

_Qadj, QadjTot, P_Swap, POver, ΣPOver, MaxEnergy = water_regulation(Nominations, Qref, T)
Qadj = Dict{Reservoir, Float64}(r => 0.0 for r in R)
for r in R
    for res in keys(_Qadj)
        if r.dischargepoint == res.dischargepoint
            Qadj[r] = _Qadj[res]
        end
    end
end
QnomO = OtherNominations(J, R,  Nominations, Qadj)

# --------------------- ShortTermScheduling and Renomination -------------------------

ReNominations = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}((participant = j, reservoir = r) => 0.0 for r in R for j in J)
for j in J
    Qnom2 = ShortTermScheduling(R, j, y[j], sample, QnomO, Omega, P)
    for r in R
        ReNominations[(participant = j, reservoir = r)] = Qnom2[r]
    end
end


_Qadj2, QadjTot2, P_Swap2, POver2, ΣPOver2, MaxEnergy2 = water_regulation(ReNominations, Qref, T)

Qadj2 = Dict{Reservoir, Float64}(r => 0.0 for r in R)
for r in R
    for res in keys(_Qadj2)
        if r.dischargepoint == res.dischargepoint
            Qadj2[r] = _Qadj2[res]
        end
    end
end
# ---------------------------- Real Time Balancing -------------------------------------

Upregulation = Dict{Participant, Vector{Float64}}(j => [0.0] for j in J)
Downregulation = Dict{Participant, Vector{Float64}}(j => [0.0] for j in J)

for j in J
    model_balancing, z_up_local, z_down_local = RealTimeBalancing(R, j, Qadj2, P_Swap2[j], y[j])
    Upregulation[j] = z_up_local
    Downregulation[j] = z_down_local
end

# ---------------------------- Evaluate Solutions --------------------------------------

Revenues = CalculateRevenues(J, y, sample, Upregulation, Downregulation)
for j in J
    println("Revenues of producer $(j.name) : ", Revenues[j], " Dollar")
end