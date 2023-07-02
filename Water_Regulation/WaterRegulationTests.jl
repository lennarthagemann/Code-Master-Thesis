try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation
using Test
# Write Unit Tests for Water Regulation Procedure

filepath_systemA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/TestSystemA.json"
filepath_systemB = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/TestSystemB.json"
filepath_systemC = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/TestSystemC.json"
filepath_systemD = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/TestSystemD.json"

filepath_nomA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestNominations/SystemA/NominationA.json"
filepath_nomB = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestNominations/SystemB/NominationB.json"
filepath_nomC = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestNominations/SystemC/NominationC.json"
filepath_nomD = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestNominations/SystemD/NominationD.json"

resA, plantsA, partsA = read_data(filepath_systemA)
resB, plantsB, partsB = read_data(filepath_systemB)
resC, plantsC, partsC = read_data(filepath_systemC)
resD, plantsD, partsD = read_data(filepath_systemD)

QnomSystemA = read_nomination(filepath_nomA, partsA, resA)
QnomSystemB = read_nomination(filepath_nomB, partsB, resB)
QnomSystemC = read_nomination(filepath_nomC, partsC, resC)
QnomSystemD = read_nomination(filepath_nomD, partsD, resD)

# Qadj_A, P_Swap_A = water_regulation(QnomSystemA, Dict{Reservoir, Float64}(r => 150 for r in resA))
# Qadj_B, P_Swap_B = water_regulation(QnomSystemB, Dict{Reservoir, Float64}(r => 150 for r in resB))
# Qadj_C, P_Swap_C = water_regulation(QnomSystemC, Dict{Reservoir, Float64}(r => 150 for r in resC))
# Qadj_D, P_Swap_D = water_regulation(QnomSystemD, Dict{Reservoir, Float64}(r => 150 for r in resD))

# @test water_regulation(resA, plantsA, partsA) == 0
# @test water_regulation(resB, plantsB, partsB) == 0
# @test water_regulation(resC, plantsC, partsC) == 0
# @test water_regulation(resD, plantsD, partsD) == 0

function create_nominations(reservoirs::Array{Reservoir}, parts::Array{Participant}, n::Int64)::Array{Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}}}
    pardict = Dict(p.name => p for p in parts)
    resdict = Dict(res.dischargepoint => res for res in reservoirs)
    nominations = Array{Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}}(undef, n)
    for i in 1:n
        nominations[i] = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}(
            (participant = pardict["A"], reservoir = resdict["r1"]) => 50 + i*10,
            (participant = pardict["B"], reservoir = resdict["r1"]) => 150,
            (participant = pardict["B"], reservoir = resdict["r2"]) => 150)
    end
    return nominations
end

Nominations_A = create_nominations(resA, partsA, 10)
# for nom in Nominations_A
#     println(nom)
#     print(water_regulation(nom, Dict{Reservoir, Float64}(r => 150 for r in resA))[2])
# end

pardictA = Dict(p.name => p for p in partsA)
resdictA = Dict(res.dischargepoint => res for res in resA)

pardictA["A"].balance[resdictA["r1"]] = 0

filepath_OvernomA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestNominations/SystemA/OverNominationA.json"
function simulate_nomination(Qnom, Qref)
    Qadj, QadjTot, P_Swap = water_regulation(Qnom, Qref)
    println("_________________________________________________________")
    println("Nomination: \n", Qnom)
    println("QadjA: \n", Qadj)
    println("QadjTotA: \n", QadjTot)
    println("P_Swap \n", P_Swap)
    total_P_Swap = Dict(p => sum(values(P_Swap[p])) for p in keys(P_Swap))
    println("Power Swap in total: ", total_P_Swap) 
    println("Power Swap cancles out: ", isapprox(sum(values(total_P_Swap)), 0, atol= 1e-6))
end

Overnomination_A = read_nomination(filepath_OvernomA, partsA, resA)
Nomination_C = read_nomination(filepath_nomC, partsC, resC)
simulate_nomination(Nomination_C, Dict{Reservoir, Float64}(r => 150 for r in resC))

Nomination_D = read_nomination(filepath_nomD, partsD, resD)
simulate_nomination(Nomination_D, Dict{Reservoir, Float64}(r => 150 for r in resD))
# QadjA, QadjTotA, P_SwapA = water_regulation(Overnomination_A, Dict{Reservoir, Float64}(r => 150 for r in resA))
# println("_________________________________________________________")
# println("Overnomination_A: ", Overnomination_A)
# println("QadjA: ", QadjA)
# println("QadjTotA: ", QadjTotA)
# println("P_SwapA: ", P_SwapA)