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
scenario_count = 1
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


"""
Private Helper Function to extract Nomination and Bid Curves as Arrays from symbolic SDDP solution 
"""
function _collect_solution_Bidding(sol, R, j; T = T)
    ks = sort(collect(keys(sol[2])))
    Qnom = Dict{Reservoir, Float64}(r => 0.0 for r in R)
    for r in R
        if j.participationrate[r] > 0.0
            Qnom[r] = sol[2][filter(value -> startswith(String(value), "Q") && occursin(r.dischargepoint, String(value)), ks)[1]]
        end
    end
    BidCurves = Dict(t => [sol[2][el] for el in filter(value -> startswith(String(value), "x") && endswith(String(value), ",$(t)]"), ks)] for t in 1:T)
    return Qnom, BidCurves
end

function _collect_solution_Short(sol, R, j; T = T)::Dict{Reservoir, Float64}
    Qnom = Dict{Reservoir, Float64}(r => 0.0 for r in R)
    for r in R
        if j.participationrate[r] > 0.0
            Qnom[r] = sol.controls[:Qnom][r]
        end
    end
    println("Nomination: \n", Qnom)
    return Qnom
end

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
-                                Bidding Problems                                   -
_____________________________________________________________________________________

=#

function Nonanticipatory_Bidding(
    all_res::Array{Reservoir},
    j::Participant,
    PPoints::Array{Float64},
    Omega,
    P;
    mu_up = 1.0,
    mu_down = 0.1,
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    stopping_rule = SDDP.BoundStalling(10, 1e-4),
    T = 24,
    Stages = stages,
    optimizer = CPLEX.Optimizer,
    DualityHandler = SDDP.ContinuousConicDuality())

    K_j = j.plants
    R = collect(filter(r -> j.participationrate[r] > 0.0, all_res))
    function subproblem_builder_nonanticipatory(subproblem::Model, node::Int64)
        @variable(subproblem, 0 <= x[i = 1:I+1, t = 1:T] <= sum(k.equivalent * k.spillreference for k in K_j), SDDP.State, initial_value=0)
        @variable(subproblem, 0 <= l[r = R] <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, lind[r = R], SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, u_start[k = K_j], SDDP.State, initial_value = 0, Bin)
        @variable(subproblem, 0 <= Qnom[r = R] <= max([k.spillreference for k in filter(k -> k.reservoir in find_ds_reservoirs(r), K)]...), SDDP.State, initial_value = 0)
        @variable(subproblem, y[t=1:T] >= 0)
        @variable(subproblem, d[t=1:T, k = K_j], Bin)
        @variable(subproblem, u[t=1:T, k = K_j], Bin)
        @variable(subproblem, BALANCE_INDICATOR[r = R], Bin)
        @variable(subproblem, 0 <= w[t=1:T, k = K_j] <= k.equivalent * k.spillreference)
        @variable(subproblem, z_up[t=1:T] >= 0)
        @variable(subproblem, z_down[t=1:T] >= 0)
        @variable(subproblem, 0 <= Qeff[t=1:T, k = K_j] <= k.spillreference)
        @variable(subproblem, 0 <= Qreal[t=1:T, r = R])
        @variable(subproblem, f[r = R] >= 0)
        @variable(subproblem, s[r = R] >= 0)
        @variable(subproblem, a[r = R])
        @constraint(subproblem, increasing[i = 1:I, t=1:T], x[i,t].out <= x[i+1,t].out)
        @constraint(subproblem, balance_ind[r = R], lind[r].out == lind[r].in - T * (Qnom[r].out - Qref[r])- s[r]) 
        @constraint(subproblem, nbal1[r = R], BALANCE_INDICATOR[r] => {Qnom[r].out <= Qref[r]})
        @constraint(subproblem, nbal2[r = R], !BALANCE_INDICATOR[r] => {0 <= lind[r].in})
        @constraint(subproblem, NoSpill[k = K_j], BALANCE_INDICATOR[k.reservoir] => {sum(Qnom[r_up].out for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference})
        @constraint(subproblem, watervalue[r = R], a[r] >= j.participationrate[r] * (l[r].in - l[r].out) * 0.5)
        if node == 1
            @stageobjective(subproblem, -sum( a[r] for r in R))
            @constraint(subproblem, balance_transfer[r = R], l[r].out == l[r].in - T * Qnom[r].out - s[r]) 
        else
            @constraint(subproblem, endcond[k = K_j], u_start[k].out == u[T,k])
            @constraint(subproblem, startcond[k = K_j], u_start[k].in == u[1,k])
            @constraint(subproblem, clearing[t=1:T], y[t] == sum(1* x[i,t].in +  1* x[i+1,t].in for i in 1:I))
            @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qnom[r].in)
            @constraint(subproblem, obligation[t=1:T], y[t] == sum(w[t,k] for k in K_j) + z_up[t] - z_down[t])
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qnom[r].out + f[r] * T - s[r])
            @constraint(subproblem, active[t=1:T, k=K_j], w[t,k] <= u[t,k] * k.spillreference * k.equivalent)
            @constraint(subproblem, startup[t=1:T-1, k=K_j], d[t,k] >= u[t+1,k] - u[t,k])
            @constraint(subproblem, production[t=1:T, k=K_j], w[t,k] == Qeff[t,k] * k.equivalent)
            @constraint(subproblem, realwater[t=1:T, k=K_j], Qeff[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)))
            SDDP.parameterize(subproblem, Omega, P) do om
                # We have to make sure that depending on the market clearing price, the coefficients are set accordingly.
                # The recourse action only applies to the real delivery, determined by the uncertain price. The other restricitions become inactive, else they make the problem infeasible.
                # The constraints that are relevant are maintained in Scenario_Index for every current time step.
                for r in R
                    JuMP.fix(f[r], om.inflow, force=true)
                end
                # Define Set of active variables for each hour
                I_t = Dict(t => 0 for t in 1:T)
                for t in 1:T
                    for i in 1:I
                        if (om.price[t] >= PPoints[i]) && (om.price[t] <= PPoints[i+1])
                            I_t[t] = i
                        end
                    end
                end
                # Include only active variables in stageobjective
                @stageobjective(subproblem , sum(om.price[t] * y[t] -  mu_up * z_up[t] + mu_down * z_down[t]  - S * sum(d[t,k] for k in K_j) for t in 1:T) -  sum(a[r] for r in R))
                # Fix / Deactivate constraints by setting their coefficients to appropriate values or all zero.
                for t in 1:T
                    for i in 1:I
                        if (i == I_t[t])
                            set_normalized_coefficient(clearing[t], x[i,t].in, -((om.price[t] - PPoints[i])/(PPoints[i+1] - PPoints[i])))
                            set_normalized_coefficient(clearing[t], x[i+1,t].in, -((PPoints[i+1] - om.price[t])/(PPoints[i+1] - PPoints[i])))
                        else
                            set_normalized_coefficient(clearing[t], x[i,t].in, 0)
                            set_normalized_coefficient(clearing[t], x[i+1,t].in, 0)
                        end
                    end
                end
            end
        end
        return
    end

    model = SDDP.LinearPolicyGraph(
        subproblem_builder_nonanticipatory;
        stages = Stages,
        sense = :Max,
        upper_bound = 1e5,
        optimizer = optimizer
    )

    SDDP.train(model; iteration_limit = 100, stopping_rules = [stopping_rule], duality_handler = DualityHandler, print_level = printlevel)

    rule = SDDP.DecisionRule(model; node = 1)
    sol = SDDP.evaluate(
        rule;
        incoming_state = Dict(Symbol("l[$(r.dischargepoint)]") => r.currentvolume for r in R),
        controls_to_record = [:l, :x, :Qnom],
    )
     
    Qnom, BiddingCurves = _collect_solution_Bidding(sol, R, j)
    return Qnom, BiddingCurves
end

function Anticipatory_Bidding(
    R::Array{Reservoir},
    j::Participant,
    PPoints::Array{Float64},
    Omega,
    P;
    mu_up = 1.0,
    mu_down = 0.1,
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    stall_bound = SDDP.BoundStalling(10, 1e-4),
    T = 24,
    Stages = stages,
    optimizer = CPLEX.Optimizer,
    DualityHandler = SDDP.ContinuousConicDuality())

    K_j = j.plants
    O, K_O = OtherParticipant(J, j, R)

    function subproblem_builder_anticipatory(subproblem::Model, node::Int64)
        @variable(subproblem, 0 <= x[i = 1:I+1, t = 1:T] <= sum(k.equivalent * k.spillreference for k in K_j), SDDP.State, initial_value=0)
        @variable(subproblem, 0 <= l[r = R] <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, lind[r = R], SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, u_start[k = K_j], SDDP.State, initial_value = 0, Bin)
        @variable(subproblem, 0 <= Qnom[r = R] <= max([k.spillreference for k in filter(k -> k.reservoir in find_ds_reservoirs(r), K)]...), SDDP.State, initial_value = 0)
        @variable(subproblem, y[t=1:T] >= 0)
        @variable(subproblem, d[t=1:T, k = K_j], Bin)
        @variable(subproblem, u[t=1:T, k = K_j], Bin)
        @variable(subproblem, BALANCE_INDICATOR[r = R], Bin)
        @variable(subproblem, 0 <= w[t=1:T, k = K_j] <= k.equivalent * k.spillreference)
        @variable(subproblem, z_up[t=1:T] >= 0)
        @variable(subproblem, z_down[t=1:T] >= 0)
        @variable(subproblem, 0 <= Qeff[t=1:T, k = K_j] <= k.spillreference)
        @variable(subproblem, 0 <= Qreal[t=1:T, r = R])
        @variable(subproblem, 0 <= Qadj[r = R])
        @variable(subproblem, Pswap[r = R])
        @variable(subproblem, Pover[k = K_O] >= 0)
        @variable(subproblem, f[r = R] >= 0)
        @variable(subproblem, s[r = R] >= 0)
        @variable(subproblem, a[r = R])
    
        @constraint(subproblem, increasing[i = 1:I, t=1:T], x[i,t].out <= x[i+1,t].out)
        @constraint(subproblem, balance_ind[r = R], lind[r].out == lind[r].in - T * (Qnom[r].out - Qref[r])- s[r]) 
        @constraint(subproblem, nbal1[r = R], BALANCE_INDICATOR[r] => {Qnom[r].out <= Qref[r]})
        @constraint(subproblem, nbal2[r = R], !BALANCE_INDICATOR[r] => {0 <= lind[r].in})
        @constraint(subproblem, NoSpill[k = K_j], BALANCE_INDICATOR[k.reservoir] => {sum(Qnom[r_up].out for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference})
        @constraint(subproblem, watervalue[r = R], a[r] >= j.participationrate[r] * (l[r].in - l[r].out) * 0.5)
        if node == 1
            @stageobjective(subproblem, -sum(a[r] for r in R))
            @constraint(subproblem, balance_transfer[r = R], l[r].out == l[r].in - T * Qnom[r].out - s[r]) 
        else
            @constraint(subproblem, adjustedflow[r = R], (j.participationrate[r] + O.participationrate[r]) * Qadj[r] - Qnom[r].in * j.participationrate[r] ==  O.participationrate[r])
            @constraint(subproblem, powerswap[r = R], Pswap[r] == j.participationrate[r] * (Qnom[r].in - Qadj[r]) - sum(Pover[k] for k in K_O))
            @constraint(subproblem, overnomination[k = K_O], Pover[k] >= k.equivalent * (Qadj[k.reservoir] - k.spillreference))
            @constraint(subproblem, endcond[k = K_j], u_start[k].out == u[T,k])
            @constraint(subproblem, startcond[k = K_j], u_start[k].in == u[1,k])
            @constraint(subproblem, clearing[t=1:T], y[t] == sum(1* x[i,t].in +  1* x[i+1,t].in for i in 1:I))
            @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qadj[r])
            @constraint(subproblem, obligation[t=1:T], y[t] - sum(Pswap[r] for r in R) == sum(w[t,k] for k in K_j) + z_up[t] - z_down[t])
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qadj[r] + f[r] * T - s[r])
            @constraint(subproblem, active[t=1:T, k=K_j], w[t,k] == u[t,k] * k.spillreference * k.equivalent)
            @constraint(subproblem, startup[t=1:T-1, k=K_j], d[t,k] >= u[t+1,k] - u[t,k])
            @constraint(subproblem, production[t=1:T, k=K_j], w[t,k] == Qeff[t,k] * k.equivalent)
            @constraint(subproblem, realwater[t=1:T, k=K_j], Qeff[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)))
            SDDP.parameterize(subproblem, Omega, P) do om
                # We have to make sure that depending on the market clearing price, the coefficients are set accordingly.
                # The recourse action only applies to the real delivery, determined by the uncertain price. The other restricitions become inactive, else they make the problem infeasible.
                # The constraints that are relevant are maiintained in Scenario_Index for every current time step.
                for r in R
                    JuMP.fix(f[r], om.inflow, force=true)
                    JuMP.set_normalized_rhs(adjustedflow[r], O.participationrate[r] * om.nomination[r])
                end
                # Define Set of active variables for each hour
                I_t = Dict(t => 0 for t in 1:T)
                for t in 1:T
                    for i in 1:I
                        if (om.price[t] >= PPoints[i]) && (om.price[t] <= PPoints[i+1])
                            I_t[t] = i
                        end
                    end
                end
                # Include only active variables in stageobjective
                @stageobjective(subproblem ,sum(om.price[t] * y[t] -  mu_up * z_up[t] + mu_down * z_down[t]  - S * sum(d[t,k] for k in K_j) for t in 1:T) - sum(a[r] for r in R))
                # Fix / Deactivate constraints by setting their coefficients to appropriate values or all zero.
                for t in 1:T
                    for i in 1:I
                        if (i == I_t[t])
                            set_normalized_coefficient(clearing[t], x[i,t].in, -((om.price[t] - PPoints[i])/(PPoints[i+1] - PPoints[i])))
                            set_normalized_coefficient(clearing[t], x[i+1,t].in, -((PPoints[i+1] - om.price[t])/(PPoints[i+1] - PPoints[i])))
                        else
                            set_normalized_coefficient(clearing[t], x[i,t].in, 0)
                            set_normalized_coefficient(clearing[t], x[i+1,t].in, 0)
                        end
                    end
                end
            end
        end
        return
    end

    model_ant = SDDP.LinearPolicyGraph(
        subproblem_builder_anticipatory;
        stages = Stages,
        sense = :Max,
        upper_bound = 1e5,
        optimizer = CPLEX.Optimizer
    )

    SDDP.train(model_ant; iteration_limit=100, stopping_rules = [stopping_rule], duality_handler = DualityHandler)

    rule_ant = SDDP.DecisionRule(model_ant; node = 1)
    sol_ant = SDDP.evaluate(
        rule_ant;
        incoming_state = Dict(Symbol("l[$(r.dischargepoint)]") => r.currentvolume for r in R),
        controls_to_record = [:l, :x, :Qnom],
    )
    
    Qnom, BiddingCurves = _collect_solution_Bidding(sol_ant, R, j)
    return Qnom, BiddingCurves
end



#=
_____________________________________________________________________________________
-                    Short-Term Optimization & Renominations                        -
_____________________________________________________________________________________

=#


function ShortTermScheduling(
    R::Array{Reservoir},
    j::Participant,
    y::Dict{Int64, Float64},
    price::Vector{Float64},
    QnomO::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64},
    Omega,
    P;
    mu_up = 1.0,
    mu_down = 0,
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    stopping_rule = [SDDP.BoundStalling(10, 1e-4)],
    T = 24,
    Stages = stages,
    optimizer = CPLEX.Optimizer)

    K_j = j.plants
    O, K_O = OtherParticipant(J, j, R)

    function subproblem_builder_short(subproblem::Model, node::Int64)
        # State Variables
        @variable(subproblem, 0 <= l[r = R] <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, lind[r = R], SDDP.State, initial_value = j.individualreservoir[r])
        @variable(subproblem, u_start[k = K_j], SDDP.State, initial_value = 1, Bin)
        # Control Variables
        @variable(subproblem, 0 <= Qnom[r = R])# <= max([k.spillreference for k in filter(k -> k.reservoir in find_ds_reservoirs(r), K)]...))
        @variable(subproblem, Qadj[r = R] >= 0)
        @variable(subproblem, d[t = 1:T, k = K_j], Bin)
        @variable(subproblem, u[t = 1:T, k = K_j], Bin)
        @variable(subproblem, BALANCE_INDICATOR[r = R], Bin)
        @variable(subproblem, 0 <= w[t = 1:T, k = K_j] <= k.equivalent * k.spillreference)
        @variable(subproblem, 0 <= Qeff[t = 1:T, k = K_j] <= k.spillreference)
        @variable(subproblem, 0 <= Qreal[t = 1:T, r = R])
        @variable(subproblem, Pswap[r = R])
        @variable(subproblem, Pover[k = K_O] >= 0)
        @variable(subproblem, s[r = R] >= 0)
        @variable(subproblem, a[r = R])
        # Random Variables
        @variable(subproblem, f[r = R] >= 0)
        # Transition Function
        if node == 1
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qadj[r] - s[r])
        else
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qadj[r] + f[r] * T - s[r])
        end
        @constraint(subproblem, balance_ind[r = R], lind[r].out == lind[r].in - T * (Qnom[r] - Qref[r])- s[r]) 
        @constraint(subproblem, startcond[k = K_j], u_start[k].in == u[1,k])
        @constraint(subproblem, endcond[k = K_j], u_start[k].out == u[T,k])
        # Constraints
        if node == 1
            @variable(subproblem, z_up[t = 1:T] >= 0)
            @variable(subproblem, z_down[t = 1:T] >= 0)
            @constraint(subproblem, obligation[t = 1:T], y[t]  == sum(w[t,k] for k in K_j) + sum(Pswap[r] for r in R) + z_up[t] - z_down[t])
        end
        @constraint(subproblem, nbal1[r = R], BALANCE_INDICATOR[r] => {Qnom[r] <= Qref[r]})
        @constraint(subproblem, nbal2[r = R], !BALANCE_INDICATOR[r] => {0 <= lind[r].in})
        @constraint(subproblem, NoSpill[r = R], BALANCE_INDICATOR[r] => {sum(Qnom[r_up] for r_up in find_us_reservoir(r)) <= min([k.spillreference for k in filter(k -> k.reservoir == r, K)]...)})
        @constraint(subproblem, adjustedflow[r = R], Qadj[r] == (Qnom[r] * j.participationrate[r] + QnomO[(participant = j, reservoir = r)] * O.participationrate[r]) / (j.participationrate[r] + O.participationrate[r]))
        @constraint(subproblem, powerswap[r = R], Pswap[r] == j.participationrate[r] * (Qnom[r] - Qadj[r]) - sum(Pover[k] for k in K_O))
        @constraint(subproblem, overnomination[k = K_O], Pover[k] >= k.equivalent * (Qadj[k.reservoir] - k.spillreference))
        @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qadj[r])
        @constraint(subproblem, active[t = 1:T, k = K_j], w[t,k] <= u[t,k] * k.spillreference * k.equivalent)
        @constraint(subproblem, startup[t = 1:T-1, k = K_j], d[t,k] >= u[t+1,k] - u[t,k])
        @constraint(subproblem, production[t = 1:T, k = K_j], w[t,k] <= Qeff[t,k] * k.equivalent)
        @constraint(subproblem, realwater[t = 1:T, k = K_j], Qeff[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)))
        @constraint(subproblem, watervalue[r = R], a[r] >= j.participationrate[r] * (lind[r].in - lind[r].out) * 0.5)
        if node > 1
            SDDP.parameterize(subproblem, Omega, P) do om
                # We have to make sure that depending on the market clearing price, the coefficients are set accordingly.
                # The recourse action only applies to the real delivery, determined by the uncertain price. The other restricitions become inactive, else they make the problem infeasible.
                # The constraints that are relevant are maiintained in Scenario_Index for every current time step.
                for r in R
                    JuMP.fix(f[r], om.inflow, force=true)
                end
                # Include only active variables in stageobjective
                if node == 1
                    # Fixed Price and Obligation
                    @stageobjective(subproblem, sum(price[t] * y[t] + mu_down * z_down[t] - mu_up * z_up[t] - S * sum(d[t,k] for k in K_j) for t in 1:T) - sum(a[r] for r in R))
                else
                    @stageobjective(subproblem, sum(om.price[t] * (sum(w[t,k] for k in K_j) + sum(Pswap[r] for r in R)) - S * sum(d[t,k] for k in K_j) for t in 1:T) - sum(a[r] for r in R))
                end
            end
        end
    return
    end

    model_short = SDDP.LinearPolicyGraph(
        subproblem_builder_short,
        stages = Stages,
        sense = :Max,
        upper_bound = sum(sum(sum(k.spillreference * k.equivalent * mu_up for k in K_j) for t in 1:T) for stage in 1:Stages),
        optimizer  = optimizer
    )   

    SDDP.train(model_short; iteration_limit = 10, stopping_rules = stopping_rule)

    rule_short = SDDP.DecisionRule(model_short; node = 1)

    sol_short = SDDP.evaluate(
        rule_short;
        incoming_state = Dict(Symbol("l[$(r.dischargepoint)]") => r.currentvolume for r in R),
        controls_to_record = [:Qnom],
    )
    Qnom = _collect_solution_Short(sol_short, R, j)
    return Qnom
end

#=
_____________________________________________________________________________________
-                              Real Timer Balancing                                 -
_____________________________________________________________________________________

=#

"""
Deterministic model to solve final optimization problem where an optimal usage of power plants and balancing market is looked for.
Returns the balancing actions, which are then used to calculate the total profits made over the planning day.
"""
function RealTimeBalancing(
    R::Vector{Reservoir},
    j::Participant,
    Qadj::Dict{Reservoir, Float64},
    Pswap::Dict{Reservoir, Float64},
    y::Dict{Int64, Float64};
    S = 0.2,
    T = 24,
    mu_up = 1.0,
    mu_down = 0.01)
    K_j = j.plants
    model_balancing = JuMP.Model(CPLEX.Optimizer)
    # Variables
    @variable(model_balancing, 0 <= z_up[t = 1:T])
    @variable(model_balancing, 0 <= z_down[t = 1:T])
    @variable(model_balancing, 0 <= w[t = 1:T, k = K_j] <= k.equivalent * k.spillreference)
    @variable(model_balancing, u[t = 1:T, k = K_j], Bin)
    @variable(model_balancing, d[t = 1:T, k = K_j], Bin)
    @variable(model_balancing, Qreal[t = 1:T, r = R])

    # Constraints
    @constraint(model_balancing, obligation[t = 1:T], y[t] == sum(w[t,k] for k in K_j) + sum(Pswap[r] for r in R) + z_up[t] - z_down[t])
    @constraint(model_balancing, adjustedflow[r = R], sum(Qreal[t, r] for t in 1:T) == T * Qadj[r])
    @constraint(model_balancing, production[t = 1:T, k = K_j], w[t,k] <= sum(Qreal[t, r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
    @constraint(model_balancing, unit[t = 1:T, k = K_j], w[t,k] <= u[t,k] * k.equivalent * k.spillreference)
    @constraint(model_balancing, startup[t = 2:T, k = K_j], d[t ,k] >= u[t,k] - u[t-1,k])
    # Objective
    @objective(model_balancing, Min, sum(mu_up * z_up[t] - mu_down * z_down[t] + sum(S * d[t, k] for k in K_j) for t in 1:T))
    optimize!(model_balancing)
    println(objective_value(model_balancing))
    return model_balancing, value.(z_up), value.(z_down)
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