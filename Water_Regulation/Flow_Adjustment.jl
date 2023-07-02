using .WaterRegulation

res1dp = "Wannsee"
res1TV = 6 * 10^12
res1CV = 3 * 10^12
res1MV = res1TV * 0.9

res1 = Reservoir(res1dp, res1TV, res1CV, res1MV)

res2dp = "Königssee"
res2TV = 5000000
res2CV = 20000
res1MV = res2TV * 0.9 

res2 = Reservoir(res2dp, res2TV, res2CV, res1MV)

hp1name = "Superkraftwerk"
hp1res= res1
hp1eq = 0.4
hp1spill = 150 

hp1 = HydropowerPlant(hp1name, hp1res, hp1eq, hp1spill)

par1name = "RWE"
par1part::Dict{Reservoir, Float16} = Dict(res1 => hp1eq)
par1plants = [hp1]
par1balance::Dict{Reservoir, Float64} = Dict(res1 => 0)

par1 = Participant(par1name, par1plants)

hp2name = "Standardkraftwerk"
hp2res= res1
hp2eq = 0.2
hp2spill = 170 

hp2 = HydropowerPlant(hp2name, hp2res, hp2eq, hp2spill)

par2name = "EON"
par2part::Dict{Reservoir, Float16} = Dict(res1 => hp2eq)
par2plants = [hp2]
par2balance::Dict{Reservoir, Float64} = Dict(res1 => 0)

par2 = Participant(par2name, par2plants)
Qnom = [Nomination(res1, par1, 120), Nomination(res1, par2, 134)]
Qref = Dict(res1 => 125)
Qadj = adjust_flow(Qnom)