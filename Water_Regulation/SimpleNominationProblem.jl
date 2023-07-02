using .WaterRegulation
using JSON
using JuMP
using CPLEX
using Distributions
using LinearAlgebra

filepath_system = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/SimpleReservoirSystem.json"
res, plants, parts = read_data(filepath_system)

j = parts[1]
K = j.plants


AUCTION_LENGTH = 24
RES_UPPER_BOUND = 1000000
Qref = 100
Qinflow = 0
Res_ref = res[1].currentvolume


model = Model(CPLEX.Optimizer)
price_lb = 0.2
price_ub = 0.5
lb_deviation = 0.9
ub_deviation = 1.1
c_Day_Ahead  = rand(Uniform(price_lb, price_ub), AUCTION_LENGTH)

@variable(model, Qnom >= 0)
@variable(model, Qeff[1:AUCTION_LENGTH, K] >= 0)
@variable(model, Qreal[1:AUCTION_LENGTH] >= 0)
@variable(model, 0 <= Res_real <= RES_UPPER_BOUND)
@variable(model, Bal)

@constraint(model,[i = 1:AUCTION_LENGTH], Qeff[i] <= Qreal[i])
for k in K
    @constraint(model, Qeff <= k.spill_reference_level)
end
@constraint(model, (sum(Qreal[i] for i in [1:AUCTION_LENGTH])/AUCTION_LENGTH == Qnom))
@constraint(model, [i = 1:AUCTION_LENGTH], lb_deviation * Qnom <= Qreal[i] <=ub_deviation * Qnom)
@constraint(model, Bal == j.balance[res[1]] + 24* (Qref - Qnom))
@constraint(model, Res_real == res[1].currentvolume - 24*(Qnom - Qinflow))
@constraint(model, Res_ref + Bal >= 0)


@objective(model, Max, (sum(dot(c_Day_Ahead, Qeff[k] * k.equivalent) for k in K))+ 0.1 * (Res_ref + Bal))

optimize!(model)