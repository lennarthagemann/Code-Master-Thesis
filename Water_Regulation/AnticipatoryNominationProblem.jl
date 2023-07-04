using .WaterRegulation
using JSON
using JuMP
using CPLEX

filepath_system = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/SimpleReservoirSystem.json"
res, plants, parts = read_data(filepath_system)

j = parts[1]
K = j.plants
partj = j.participationrate[resr[1]]

O = OtherProducer(parts, j)

QnomO = 70
partO = O.participationrate[res[1]]
"""
Solve the optimization problem of O to determine the necessary parameters for the anticipative optimization
"""
AUCTION_STEPS = 1
AUCTION_LENGTH = 24
RES_UPPER_BOUND = 1000000
Qref = 70
Qinflow = 0
Res_ref = res[1].currentvolume
M = RES_UPPER_BOUND

model = Model(CPLEX.Optimizer)
c_Day_Ahead  = fill(0.3, AUCTION_LENGTH ,1)

@variable(model, Qnomj)
@variable(model, Qeff >= 0)
@variable(model, Qadj >= 0)
@variable(mode, Qmax )
@variable(model, 0 <= Res_real <= RES_UPPER_BOUND)
@variable(model, Bal)
@variable(model, y, Bin)

@constraint(model, Qadj == (partj * Qnom  + partO * QnomO)/(partj + partO))
@constraint(model, Qnom <= (Qmax(partj + partO) - QnomO * partO)/(partj) + M*y)
@constraint(model, 0 <= (Res_ref + Bal) + M*y)
@constraint(model, Qeff <= Qadj )
for k in K
    @constraint(model, Qeff <= k.spill_reference_level)
end
@constraint(model, Bal == j.balance[res[1]] + 24* (Qref - Qnom))
@constraint(model, Res_real == res[1].currentvolume - 24*(Qnom + Qinflow))

@objective(model, Max, (sum(sum(Qeff * k.equivalent for k in K) * c_Day_Ahead)) + 0.1 * (Res_ref + Bal))

optimize!(model)