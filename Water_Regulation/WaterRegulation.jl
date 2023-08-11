module WaterRegulation

using JSON
using Printf
using JuMP
using CPLEX
using SDDP
export HydropowerPlant, Reservoir, Participant, adjust_flow!, calculate_balance, update_reservoir, update_ind_reservoir, Calculate_Ersmax, Calculate_POver, power_swap, find_us_reservoir, find_ds_reservoirs, connect_reservoirs, read_nomination, read_data, water_regulation, OtherParticipant, CalculateQmax, Calculate_Qover, partAvg, SimplePartAvg, SumPartAvg, calculate_produced_power, total_power, Nonanticipatory_Bidding, Anticipatory_Bidding, ShortTermScheduling, RealTimeBalancing
 
mutable struct Reservoir
    dischargepoint::String
    totalvolume::Float64
    currentvolume::Float64
    maxvolume::Float64
    upstreamreservoir::Union{Array{Reservoir}, Nothing}
    downstreamreservoir::Union{Reservoir, Nothing}
    function Reservoir(dischargepoint::String, totalvolume::Float64, currentvolume::Float64, maxvolume::Float64)
        if totalvolume < maxvolume
            throw(DomainError("Total Volume of reservoir cannot be smaller than maximum dischargable volume"))
        end
        return new(dischargepoint, totalvolume, currentvolume, maxvolume, nothing, nothing)
    end
    # Write dictionary to keep track of all constructed reservoirs with their dischargepoint mapping on the reservoir
    # This is updated once a new Reservoir is constructed.
    function Reservoir()
        return new("TestReservoir", 100, 90, 90, nothing, nothing)
    end
end

"""
Connects upstream reservoirs (Array r2) with the downstream reservoir r1.
Update structures in respective reservoirs
"""
function connect_reservoirs(dr::Reservoir, ur::Array{Reservoir})
    for r in ur
        if typeof(dr.upstreamreservoir) == Nothing
            dr.upstreamreservoir = [r]
        else
            push!(dr.upstreamreservoir, r)
        end
        r.downstreamreservoir = dr
    end
end

"""
Find all upstream reservoirs of a given discharge point, including itself. Relevant for the calculation of any participation rate,
as participants can nominate for discharge at any upstream reservoir of their powerplant.
"""
function find_us_reservoir(r::Reservoir)::Array{Reservoir}
    queue::Array{Reservoir} = [r]
    usres::Array{Reservoir} = []
    # Write while loop that terminates when queue is empty. Add upstream_reservoirs of elements in queue during iteration.
    while !isempty(queue)
        currentres::Reservoir = popfirst!(queue)
        @assert (typeof(currentres) == Reservoir)
        push!(usres, currentres)
        if typeof(currentres.upstreamreservoir) != Nothing
            for res in currentres.upstreamreservoir
                push!(queue, res)
            end 
        end
    end
    return usres
end

function find_ds_reservoirs(r::Reservoir)::Array{Reservoir}
    queue::Array{Reservoir} = [r]
    dsres::Array{Reservoir} = []
    while !isempty(queue)
        currentres::Reservoir = popfirst!(queue)
        @assert (typeof(currentres) == Reservoir)
        push!(dsres, currentres)
        if typeof(currentres.downstreamreservoir) != Nothing
            push!(queue, currentres.downstreamreservoir)
        end
    end
    return dsres
end

struct HydropowerPlant
    name::String
    reservoir::Reservoir
    equivalent::Float64
    spillreference::Float64
    function HydropowerPlant(name::String, reservoir::Reservoir, equivalent::Float64, spillreference::Float64)
        return new(name, reservoir, equivalent, spillreference)
    end
    function HydropowerPlant()
        return new("TestPlant", Reservoir(), 0.1, 1)
    end
end

struct Participant
    name::String
    plants::Array{HydropowerPlant}
    participationrate::Dict{Reservoir, Float64} 
    individualreservoir::Dict{Reservoir, Float64}
    function calculate_participation(plants, res)
        prate = Dict{Reservoir, Float64}(r => 0.0 for r in res)
        for p in plants
            for res in find_us_reservoir(p.reservoir)
                prate[res] += p.equivalent
            end  
        end
        return prate::Dict{Reservoir, Float64}
    end
    function Participant(name::String, plants::Array{HydropowerPlant}, res::Array{Reservoir})
        prate = calculate_participation(plants, res)
        ind_res  = Dict{Reservoir, Float64}(r => r.currentvolume for r in res)
        return new(name, plants, prate, ind_res)   
    end
    function Participant(name::String, plants::Array{HydropowerPlant}, ind_res::Dict{Reservoir, Float64}, res::Array{Reservoir})
        prate = calculate_participation(plants, res)
        return new(name, plants, prate, ind_res)   
    end
    function Participant()
        return new("Default", [], Dict{Reservoir, Float64}(), Dict{Reservoir, Float64}()) 
    end
end

"""
Return an aggregated other producer: From the point of one producer p, return a new producer who owns the complement of all power in the river system.
"""
function OtherParticipant(parts::Array{Participant}, p::Participant, res::Array{Reservoir})
    O = filter(x -> x != p, parts)
    plantsO::Array{HydropowerPlant} = vcat([part.plants for part in O]...)
    pO = Participant("Other", plantsO, res)
    return pO, plantsO
end 


function Base.show(io::IO, hp::HydropowerPlant)
    print(io, hp.name)
end

function Base.show(io::IO, r::Reservoir)
    print(io, r.dischargepoint)
end

function Base.show(io::IO, p::Participant)
    print(io, p.name)
end

function Base.show(io::IO, Qadj_All::Dict{Reservoir, Float64})
    println(io, "________________________________")
    for (k,v) in Qadj_All
        println(io, "$(rpad(k.dischargepoint, 8)) | $(rpad(v, 8))")
    end
end

function Base.show(io::IO, Qnom::Dict{Tuple{Participant,Reservoir}, Float64})
    println(io, "Nominations of all producers at every discharge point")
    for (k,v) in Qnom
        println(io, "________________________________")
        for el in v
            println(io, el)
        end
    end
end

function Base.show(io::IO, POver::Dict{Participant, Dict{HydropowerPlant, Float64}})
    for (k,v) in POver
        println(k.name)
        println("______________________________")
        for (k1,v1) in v
            println(io, "$(rpad(k1.name, 8)) | $(rpad(v1, 8))")
        end
    end
end

"""
Calculates the balance of all producer at a discharge point after nomination.
Sets the new balance on the participant objects
"""
function calculate_balance(Qref::Dict{Reservoir, Float64}, Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, d::Reservoir)
    for (nom,value) in Qnom
        nom.part.individualreservoir[nom.res] += (Qref[nom.res] - value)*24
    end
end

"""
Calculate the adjusted flow at every discharge point. The flow aims to preserve the energy generated from individual nominations.
We iteratively add all nominations, multiplied with their participation rate. The sum is then divided by the sum of all participation rates.
If a producer has a negative balance, he has to reduce his bid so that the adjusted flow is less than Qmax.
"""
function adjust_flow!(
    Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64},
    Qmax::Dict{HydropowerPlant, Float64}, 
    adjust_nomination=false
    )
    Qadj::Dict{Reservoir, Float64} = Dict{Reservoir, Float64}(r => 0 for r in unique([nom.reservoir for nom in keys(Qnom)]))
    # Qnom_cleaned = adjusted_nominations(Qnom)
    prates_res = Dict{Reservoir, Float64}(r => 0 for r in unique([nom.reservoir for nom in keys(Qnom)]))
    for (nom,v) in Qnom
        @assert haskey(nom.participant.participationrate, nom.reservoir)
        Qadj[nom.reservoir] += v * nom.participant.participationrate[nom.reservoir]
        prates_res[nom.reservoir] += nom.participant.participationrate[nom.reservoir]
    end
    for r in keys(Qadj)
        Qadj[r] = Qadj[r] / prates_res[r]
    end
    if adjust_nomination == true # Adjust a nomination if a producer has a negativebalance
        Qmaxr = Dict{Reservoir, Float64}(r => min([Qmax[k] for k in filter(x -> x.reservoir == r, keys(Qmax))]...) for r in unique([nom.reservoir for nom in keys(Qnom)]))
        for r in unique([nom.reservoir for nom in keys(Qnom)])
            partSum = sum(n.participant.participationrate[r] for n in filter(x -> x.reservoir == r ,keys(Qnom)))
            for (nom,value) in filter(kv -> kv[1].reservoir == r, Qnom)
                QnomO = sum(val for (key,val) in filter(x -> x[1].reservoir == r ,Qnom)) - value
                if Qadj[r] > Qmaxr[r] && !(nom.participant.individualreservoir[nom.reservoir]< 0) 
                    partj = nom.part.participationrate[r]
                    partO = partSum - partj
                    @info "The participant $(nom.participant.name) has a negative balance at $(r.dischargepoint) and his nomination of $(nom.value)has to be adjusted to $((Qmaxr[r] * partSum - QnomO*partO)/partj) \n"
                    value = (Qmaxr[r] * partSum - QnomO*partO)/partj
                    Qadj[r] = sum((val * nom.participant.participationrate[r])/(partSum) for (nom, val) in filter(x -> x[1].reservoir == r ,Qnom))
                    @assert Qadj[r] <= Qmaxr[r] 
                end
            end
        end
    end
    # Get every combination of producer and discharge point from the tuple keys of Qnom, including combinations that are not in the set of keys.
    part_entries = unique(collect(map(x -> x[1], collect(keys(Qnom)))))
    res_entries = unique(collect(map(x -> x[2], collect(keys(Qnom)))))
    # Nomination is equal to adjusted flow if participant has no participation rate
    for p in part_entries
        for r in res_entries
            if (!haskey(Qnom, (participant = p, reservoir = r)))
                Qnom[(participant = p, reservoir = r)] = Qadj[r]
            end
        end
    end
    # Obtain the total nominated and adjusted discharge. Add up the values of upstream reservoirs (including this one). 
    QadjTot::Dict{Reservoir, Float64} = Dict{Reservoir, Float64}(r => sum([Qadj[us] for us in find_us_reservoir(r)]) for r in keys(Qadj))
    # This is a little trickier for the nominations: Here we have to find the upstream nominations of each producer and add them up.
    # So: Find the upstream reservoirs and the nominations of each producer at these reservoirs. Add them up.
    QnomTot = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}() 
    for k in keys(Qnom)
        QnomTot[(participant = k.participant, reservoir = k.reservoir)] = sum([Qnom[(participant = k.participant, reservoir = us)] for us in find_us_reservoir(k.reservoir)])
    end 
    for p in part_entries
        for r in res_entries
            @assert haskey(Qnom, (participant = p, reservoir = r))
            @assert haskey(QnomTot, (participant = p, reservoir = r))
        end
    end
    return Qadj, Qnom, QadjTot, QnomTot
end

"""
Update the real (!) reservoir volume after adjustment of flow 
"""
function update_reservoir(Qadj::Dict{Reservoir, Float64}, T::Int64; round_decimal_places = 5)
   for (d, adj) in Qadj
        d.currentvolume -= adj*T
    end
end

"""
Update the individual reservoir by updating the balance of every participant through their new nomination.
Every participant has a field balance, it is the updated by the difference of nomination value and reference flow at every reservoir.
(Times 24, as it is done for the entire day)
"""
function update_ind_reservoir(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qref::Dict{Reservoir, Float64})
    for (nom, value) in Qnom
        if haskey(nom.participant.individualreservoir, nom.reservoir)
            nom.participant.individualreservoir[nom.reservoir] += (Qref[nom.reservoir] - value)*24
        end
    end
end

"""
Reads in Data about Power Producers (Participants), Hydropower Plants and Reservoirs from .json file.
Returns their corresponding objects
"""
function read_data(filename::String)
    d = JSON.parsefile(filename)
    # Unpack Participants, Rservoir and HydropowerPlants individiually
    res = [Reservoir(r["dischargepoint"], float(r["totalvolume"]), float(r["currentvolume"]), float(r["maxvolume"])) for r in d["Reservoirs"]]
    # Identify Reservoir by its dischargepoint
    resdict = Dict(r.dischargepoint => r for r in res)
    for (lr, ur) in d["Connections"]
        connect_reservoirs(resdict[lr], [resdict[r] for r in ur])
    end
    # Unpack HydropowerPlant
    plants = [HydropowerPlant(p["name"], resdict[p["reservoir"]], float(p["equivalent"]), float(p["spillreference"])) for p in d["HydropowerPlants"]]
    # Unpack Participants
    plantdict = Dict(p.name => p for p in plants)
    parts = [Participant(part["name"], [plantdict[pname] for pname in part["plants"]],  Dict(resdict[k] => float(v) for (k,v) in part["individualreservoir"]), res) for part in d["Participants"]]
    # Set participation rate to zero for missing reservoirs in the participants participation rate.
    for r in res
        for p in parts
            if !(haskey(p.participationrate, r))
                p.participationrate[r] = 0.0
            end
        end
    end
    return res::Array{Reservoir}, plants::Array{HydropowerPlant}, parts::Array{Participant}
end

"""
Read in nomination. After power participant solve their optimization problem, their nominations are saved as a dictionary
with their name (just their name!) as key and the nomination as value, which is contained in a dictionary for all discharge points. 
This is done because nomionations, adjustment, balance updates and power swaps are calculated for each discharge point individually.
The output should be usable in !adjust_flow.
"""
function read_nomination(filename::String, participants::Array{Participant}, reservoirs::Array{Reservoir})::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}
    pardict = Dict(p.name => p for p in participants)
    resdict = Dict(res.dischargepoint => res for res in reservoirs)
    d = JSON.parsefile(filename)
    Qnom = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}((participant = pardict[el[1]], reservoir = resdict[el[3]]) => el[2] for el in d)
    return Qnom
end
"""
Weighted Average participation rate of a producer upstream (and including) of a given reservoir: Using this participation rate
we can calculate the sum of adjusted flows as one adjusted flow with the total nomination.
"""
function partAvg(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, QnomTot::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64})
    part_entries = unique(collect(map(x -> x[1], collect(keys(Qnom)))))
    res_entries = unique(collect(map(x -> x[2], collect(keys(Qnom)))))
    pAvg = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}()
    ΣpAvg = Dict{Reservoir, Float64}()
    for r in res_entries
        us_res = find_us_reservoir(r)
        for p in part_entries
            pAvg[(participant = p, reservoir = r)] = sum(p.participationrate[us] * Qnom[(participant = p, reservoir = us)] for us in us_res) / QnomTot[(participant = p, reservoir = r)]
        end
        ΣpAvg[r] = sum(pAvg[(participant = p, reservoir = r)] for p in part_entries) 
    end
    return pAvg, ΣpAvg
end
"""
Sum of all average participation rates of all producers at a discharge point, used for calculation of overnomination in water
"""
function SumPartAvg(r::Reservoir, ps::Array{Participant})::Float64
    return sum(partAvg(p, r) for p in ps)
end

"""
Unfortunately I can't use the (hopefully real) pAvg calculation, as it introduces nonlinearities in the model.
Here is a simpler version, which doesn't the nomination values of the producers, but only their participation rates. 
"""
function SimplePartAvg(parts::Array{Participant}, res::Array{Reservoir})
    pAvg = Dict{Participant, Dict{Reservoir, Float64}}()
    ΣpAvg = Dict{Reservoir, Float64}()
    for p in parts
        pAvg[p] = Dict{Reservoir, Float64}(r => sum(p.participationrate[us] for us in find_us_reservoir(r))/length(find_us_reservoir(r)) for r in res)
    end
    for r in res
        ΣpAvg[r] = sum(pAvg[p][r] for p in parts)
    end
    return pAvg, ΣpAvg
end

"""
Function that calculates the maximum admissable flow at every HydropowerPlant. Necessary to calculate the reduction of power swaps.
"""
function CalculateQmax(QnomTot::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qref::Dict{Reservoir, Float64})::Dict{HydropowerPlant, Float64}
    Qmax = Dict{HydropowerPlant, Float64}()
    for (nom, value) in QnomTot
        for p in nom.participant.plants
            if p.reservoir == nom.reservoir
                Qmax[p] =  max(min(p.spillreference, Qref[p.reservoir]), value)
            end
        end
    end
    return Qmax
end

"""
Maximum compensated energy, real amount of energy lost through spillage
"""
function Calculate_Ersmax(plants::Array{HydropowerPlant}, QadjTot::Dict{Reservoir, Float64})::Dict{HydropowerPlant, Float64}
    Ersmax = Dict{HydropowerPlant, Float64}()
    for k in filter(plant -> plant in plants, plants)
        Ersmax[k] = max((QadjTot[k.reservoir] - k.spillreference) * k.equivalent, 0)
    end
    return Ersmax
end
"""
Calculate the overnomination at every discharge point in physical power. Necessary to calculate the reduction of power swaps.
The Overnomination at every plant that the producer doesn't own is summed to get the overnomination at every discharge point.
We use the partAvg and SumPartAvg to evaluate how much water is actually overnominated by a single producer.
"""
function Calculate_POver(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, QnomTot::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qmax::Dict{HydropowerPlant, Float64})
    parts = unique(collect(map(x -> x[1], collect(keys(Qnom)))))
    plants = vcat([p.plants for p in parts]...)
    res = unique(collect(map(x -> x[2], collect(keys(Qnom)))))
    POver = Dict{Participant, Dict{HydropowerPlant, Float64}}(p => Dict{HydropowerPlant, Float64}(k => 0.0 for k in plants) for p in parts)
    # pAvg, ΣpAvg = partAvg(Qnom,QnomTot)
    pAvg, ΣpAvg = SimplePartAvg(parts,res)
    for k in plants
        for part in parts
            if k in part.plants
                POver[part][k] = 0
            else
                r = k.reservoir
                @assert haskey(QnomTot, (participant = part, reservoir = r)) "QnomTot has the missing key"
                POver[part][k] = max(QnomTot[(participant = part, reservoir = r)] - Qmax[k], 0) * k.equivalent * pAvg[part][r]/ΣpAvg[r]
            end
        end
    end
    ΣPOver = Dict{HydropowerPlant, Float64}(k => sum([POver[p][k] for p in parts]) for k in plants)
    return POver, ΣPOver
end

"""
Calculation of power swap including reduction. Makes sure that the reduction is fair to both the trigger and sufferer.
"""
function power_swap(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qadj_All::Dict{Reservoir, Float64}, 
    POver::Dict{Participant, Dict{HydropowerPlant, Float64}}, ΣPOver::Dict{HydropowerPlant, Float64},
    Ersmax::Dict{HydropowerPlant, Float64}, parts::Array{Participant})::Dict{Participant, Dict{Reservoir, Float64}}
    P_Swap = Dict{Participant, Dict{Reservoir, Float64}}(p => Dict{Reservoir, Float64}(r => 0.0 for r in keys(p.participationrate)) for p in parts)
    plants = vcat([p.plants for p in parts]...)
    Red, ΣRed = PS_Reduction(Ersmax, POver, ΣPOver, parts, plants)
    for (nom, value) in Qnom
        P_Swap[nom.participant][nom.reservoir] = (value - Qadj_All[nom.reservoir]) * nom.participant.participationrate[nom.reservoir] - ΣRed[nom.participant][nom.reservoir]
    end
    return P_Swap
end

"""
The Power Swap will be reduced depending on the spillage that occurs at another plant by overnomination.
The Reduction is limited by the theoretical overnomination (POver) and the real energy lost through spillage (Ersmax).
That way producers can't be compensated for more than is actually lost, and producers can't be reduced by more than they are accountable for.
The reduction has to be accounted for on both sides! What is reduced, has to increase on the other end of the line.
"""
function PS_Reduction(Ersmax::Dict{HydropowerPlant, Float64}, POver, ΣPOver, parts::Array{Participant}, plants)
    res = unique([p.reservoir for p in plants])
    Red = Dict{Participant, Dict{HydropowerPlant, Float64}}(p => Dict{HydropowerPlant, Float64}(k => 0.0 for k in plants) for p in parts)
    for k in plants
        # Get the owner of the plant
        p_k = filter(x -> k in x.plants, parts)[1]
        for p in parts
            if Ersmax[k] >= ΣPOver[k]
                Red[p][k] += POver[p][k] 
                Red[p_k][k] -= POver[p][k] 
                
            else
                Red[p][k] += POver[p][k] * Ersmax[k]/ΣPOver[k]
                Red[p_k][k] -= POver[p][k] * Ersmax[k]/ΣPOver[k]
            end
        end
    end
    plants_at_res = Dict{Reservoir, Array{HydropowerPlant}}(r => collect(unique(filter(k -> k.reservoir == r, plants))) for r in res)
    # Sum up Reduction of Powerplants to one value for each reservoir
    ΣRed::Dict{Participant, Dict{Reservoir, Float64}} = Dict(p => Dict(r => 0 for r in res) for p in parts)
    for p in parts
        for r in res
            for plant in plants_at_res[r]
                ΣRed[p][r] += Red[p][plant]
            end
        end
    end
    return Red, ΣRed
end

"""
Calculate how much power is produced hourly at every power plant.
Take the efficient flow of water and multiply with the power plants' efficiency.
Return an array of hourly prices.
"""
function calculate_produced_power(Qeff::Dict{Reservoir, Float64}, parts::Array{Participant})::Dict{Participant, Float64}
    return Dict{Participant, Float64}(p => sum([Qeff[k.reservoir] * k.equivalent for k in p.plants]) for p in parts)
end

"""
Calculate the total amount of power generated by every participant after adjustment and power swap (as well as reduction)
The total amount of power is equivalent to the sum of power generated at every haydropowerplant and the power swap at that dischargepoint
"""
function total_power(Qadj::Dict{Reservoir, Float64}, P_Swap::Dict{Participant, Dict{Reservoir, Float64}}, parts::Array{Participant})::Dict{Participant, Float64}
    Produced_Power = calculate_produced_power(Qadj, parts)
    return Dict{Participant, Float64}(p => Produced_Power[p] + sum(values(P_Swap[p])) for p in parts)
end

"""
This is the main function: here the entire water regulation process in one step takes place.
- Nominations are taken as input
- The adjusted flow is calculated
- The power swap and reduction of it are calculated
- real and individual reservoirs are updated

"""
function water_regulation(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qref::Dict{Reservoir, Float64}, T::Int64)
    parts = unique([nom.participant for nom in keys(Qnom)])
    plants = vcat([part.plants for part in parts]...)
    QnomTot_prev = Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}() 
    for k in keys(Qnom)
        QnomTot_prev[(participant = k.participant, reservoir = k.reservoir)] = sum([Qnom[(participant = k.participant, reservoir = us)] for us in find_us_reservoir(k.reservoir)])
    end 
    Qmax = CalculateQmax(QnomTot_prev, Qref)
    Qadj, Qnom, QadjTot, QnomTot = adjust_flow!(Qnom, Qmax)
    update_reservoir(Qadj, T)
    update_ind_reservoir(Qnom, Qref)
    MaxEnergy = Calculate_Ersmax(plants, QadjTot)
    POver, ΣPOver = Calculate_POver(Qnom, QnomTot, Qmax)
    P_Swap = power_swap(Qnom, Qadj, POver, ΣPOver, MaxEnergy, parts)
    # Total_Power = total_power(Qadj_All, P_Swap)
    return Qadj, QadjTot, P_Swap, POver, ΣPOver, MaxEnergy
end



# --------------------------- Optimization Functions from SDDP -------------------------------------- 

"""
Helper Function to extract solution from Bidding Problem
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
"""
Helper Function to extract solution from ShortTermScheduling Problem.
"""
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
        @variable(subproblem, BALANCE_INDICATOR[r = R], Bin)
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
            @constraint(subproblem, start_transfer[k = K_j] , u_start[k].out == u_start[k].in)
        else
            @variable(subproblem, 0 <= w[t=1:T, k = K_j] <= k.equivalent * k.spillreference)
            @variable(subproblem, z_up[t=1:T] >= 0)
            @variable(subproblem, z_down[t=1:T] >= 0)
            @variable(subproblem, 0 <= Qreal[t=1:T, r = R])
            @variable(subproblem, f[r = R] >= 0)
            @variable(subproblem, y[t=1:T] >= 0)
            @variable(subproblem, d[t=1:T, k = K_j], Bin)
            @variable(subproblem, u[t=1:T, k = K_j], Bin)
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qnom[r].out + f[r] * T - s[r])
            @constraint(subproblem, endcond[k = K_j], u_start[k].out == u[T,k])
            @constraint(subproblem, startcond[k = K_j], u_start[k].in == u[1,k])
            @constraint(subproblem, clearing[t=1:T], y[t] == sum(1* x[i,t].in +  1* x[i+1,t].in for i in 1:I))
            @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qnom[r].in)
            @constraint(subproblem, obligation[t=1:T], y[t] == sum(w[t,k] for k in K_j) + z_up[t] - z_down[t])
            @constraint(subproblem, active[t=1:T, k=K_j], w[t,k] <= u[t,k] * k.spillreference * k.equivalent)
            @constraint(subproblem, startup[t=1:T-1, k=K_j], d[t,k] >= u[t+1,k] - u[t,k])
            @constraint(subproblem, production[t=1:T, k=K_j], w[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
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
            @constraint(subproblem, startup[t=1:T-1, k=K_j], d[t,k] >= u[t+1,k] - u[t,k])
            @constraint(subproblem, active[t=1:T, k=K_j], w[t,k] <= u[t,k] * k.spillreference * k.equivalent)
            @constraint(subproblem, production[t=1:T, k=K_j], w[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
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
        
        SDDP.train(model_ant; iteration_limit=100, stopping_rules = [stopping_rule], duality_handler = DualityHandler, print_level = printlevel)
        
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
    mu_down = 0.1,
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
        @variable(subproblem, u_start[k = K_j], SDDP.State, initial_value = 0, Bin)
        # Control Variables
        @variable(subproblem, 0 <= Qnom[r = R] <= 10000)
        @variable(subproblem, d[t = 1:T, k = K_j], Bin)
        @variable(subproblem, u[t = 1:T, k = K_j], Bin)
        @variable(subproblem, BALANCE_INDICATOR[r = R], Bin)
        @variable(subproblem, 0 <= w[t = 1:T, k = K_j] <= k.equivalent * k.spillreference)
        @variable(subproblem, 0 <= Qreal[t = 1:T, r = R])
        @variable(subproblem, s[r = R] >= 0)
        @variable(subproblem, a[r = R])
        # Random Variables
        @variable(subproblem, f[r = R] >= 0)
        @constraint(subproblem, balance_ind[r = R], lind[r].out == lind[r].in - T * (Qnom[r] - Qref[r])- s[r]) 
        @constraint(subproblem, startcond[k = K_j], u_start[k].in == u[1,k])
        @constraint(subproblem, endcond[k = K_j], u_start[k].out == u[T,k])
        # Constraints
        if node == 1
            @variable(subproblem, z_up[t = 1:T] >= 0)
            @variable(subproblem, z_down[t = 1:T] >= 0)
            @variable(subproblem, Qadj[r = R] >= 0)
            @variable(subproblem, Pswap[r = R])
            @variable(subproblem, Pover[k = K_O] >= 0)
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qadj[r] - s[r])
            @constraint(subproblem, obligation[t = 1:T], y[t]  == sum(w[t,k] for k in K_j) + sum(Pswap[r] for r in R) + z_up[t] - z_down[t])
            @constraint(subproblem, powerswap[r = R], Pswap[r] == j.participationrate[r] * (Qnom[r] - Qadj[r]) - sum(Pover[k] for k in K_O))
            @constraint(subproblem, overnomination[k = K_O], Pover[k] >= k.equivalent * (Qadj[k.reservoir] - k.spillreference))
            @constraint(subproblem, adjustedflow[r = R], Qadj[r] == (Qnom[r] * j.participationrate[r] + QnomO[(participant = j, reservoir = r)] * O.participationrate[r]) / (j.participationrate[r] + O.participationrate[r]))
            @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qadj[r])
        else
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - T * Qnom[r] + f[r] * T - s[r])
            @constraint(subproblem, nomination[r = R], sum(Qreal[t,r] for t in 1:T) == T * Qnom[r])
        end
        @constraint(subproblem, nbal1[r = R], BALANCE_INDICATOR[r] => {Qnom[r] <= Qref[r]})
        @constraint(subproblem, nbal2[r = R], !BALANCE_INDICATOR[r] => {0 <= lind[r].in})
        @constraint(subproblem, NoSpill[r = R], BALANCE_INDICATOR[r] => {sum(Qnom[r_up] for r_up in find_us_reservoir(r)) <= min([k.spillreference for k in filter(k -> k.reservoir == r, K)]...)})
        @constraint(subproblem, startup[t = 1:T-1, k = K_j], d[t,k] >= u[t+1,k] - u[t,k])
        @constraint(subproblem, active[t = 1:T, k = K_j], w[t,k] <= u[t,k] * k.spillreference * k.equivalent)
        @constraint(subproblem, production[t = 1:T, k = K_j], w[t,k] <= sum(Qreal[t,r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
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
                @stageobjective(subproblem, sum(om.price[t] * sum(w[t,k] for k in K_j) - S * sum(d[t,k] for k in K_j) for t in 1:T) - sum(a[r] for r in R))
            end
        else
            # Fixed Price and delivery in first stage, only think about minimum water value loss and productions costs
            @stageobjective(subproblem, sum(Prices[1][t] * y[t] - mu_up * z_up[t] + mu_down * z_down[t] - S * sum(d[t,k] for k in K_j) for t in 1:T) - sum(a[r] for r in R))
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
        
    SDDP.train(model_short; iteration_limit = 10, stopping_rules = stopping_rule, print_level = printlevel)
    
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
-                              Real Time Balancing                                  -
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
    set_silent(model_balancing)
    optimize!(model_balancing)
    return model_balancing, value.(z_up), value.(z_down)
end

#=
_____________________________________________________________________________________
-                              Single Owner Models                                  -
_____________________________________________________________________________________

=#

function SingleOwnerBdding(
    R::Array{Reservoir},
    K::Array{HydropowerPlant},
    PPoints::Array{Float64},
    Omega,
    P;
    mu_up = 1.0,
    mu_down = 0.1,
    S = 0.3,
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    stopping_rule = [SDDP.BoundStalling(10, 1e-4)],
    T = 24,
    Stages = stages,
    optimizer = CPLEX.Optimizer,
    DualityHandler = SDDP.ContinuousConicDuality())
    function subproblem_builder_single_bidding(subproblem::Model, node::Int64)
        # State Variables
        @variable(subproblem, 0 <= l[r = R] <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, ustart[k = K], Bin, SDDP.State, initial_value = 0)
        @variable(subproblem, 0 <= x[i = 1:I+1, t = 1:T] <= sum(k.equivalent * k.spillreference for k in K), SDDP.State, initial_value = 0)
        # Transition Function
        @constraint(subproblem, increasing[i = 1:I, t=1:T], x[i,t].out <= x[i+1,t].out)
        if node == 1
            # We only concern ourselves with bidding in the first stage.
            @stageobjective(subproblem, 0)
            @constraint(subproblem, balance_transfer[r = R], l[r].out == l[r].in) # No Inflow In first stage (or fixed value)
            @constraint(subproblem, start_transfer[k = K], ustart[k].out == ustart[k].in)
        else
            # Some Constraints and variables such as production specific components only become relevant in stage >= 2
            @variable(subproblem, y[t = 1:T] >= 0)
            @variable(subproblem, a[r = R])
            @variable(subproblem, 0 <= Qreal[t = 1:T, r = R])
            @variable(subproblem, 0 <= w[t = 1:T, k = K] <= k.equivalent * k.spillreference)
            @variable(subproblem, u[t = 1:T, k = K], Bin)
            @variable(subproblem, d[t = 1:T, k = K], Bin)
            @variable(subproblem, s[r = R] >= 0)
            @variable(subproblem, z_up[t = 1:T] >= 0)
            @variable(subproblem, z_down[t = 1:T] >= 0)
            # Random Variables
            @variable(subproblem, f[r = R])
            @constraint(subproblem, watervalue[r = R], a[r] >= sum(k.equivalent for k in filter( k-> k.reservoir in find_ds_reservoirs(r), K)) * (l[r].in - l[r].out) * 0.5)
            @constraint(subproblem, clearing[t=1:T], y[t] == sum(1* x[i,t].in +  1* x[i+1,t].in for i in 1:I-1))
            @constraint(subproblem, balance[r = R], l[r].out == l[r].in - sum(Qreal[t,r] for t in 1:T) + T * f[r] - s[r])
            @constraint(subproblem, planttrans1[k = K], ustart[k].in == u[1,k])
            @constraint(subproblem, planttrans2[k = K], ustart[k].out == u[T,k])
            # Constraints
            @constraint(subproblem, startup[t = 1:T-1, k = K], d[t,k] >= u[t+1,k] - u[t,k])
            @constraint(subproblem, activeplant[t = 1:T, k = K], w[t,k] <= u[t,k] * k.equivalent * k.spillreference)
            @constraint(subproblem, production[t = 1:T, k = K], w[t,k] <= sum(Qreal[t, r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
            @constraint(subproblem, obligation[t = 1:T], y[t] == sum(w[t, k] for k in K) + z_up[t] - z_down[t])
            # Parameterize Uncertainty and Objective Function
            SDDP.parameterize(subproblem, Omega, P) do om
                for r in R
                    JuMP.fix(f[r], om.inflow[r], force=true)
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
                # Parameterize objective through uncertain price
                @stageobjective(subproblem, sum(om.price[t] * y[t] -  mu_up * z_up[t] + mu_down * z_down[t]  - sum(S * d[t,k] for k in K) - sum(a[r] for r in R) for t in 1:T))
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
    end
    model_single_bidding = SDDP.LinearPolicyGraph(subproblem_builder_single_bidding; stages = Stages, sense = :Max,
    upper_bound = Stages * sum(sum(k.spillreference * k.equivalent for k in K) for t in 1:T for s in 1:Stages), optimizer = optimizer)
    SDDP.train(model_single_bidding; stopping_rules = stopping_rule, iteration_limit = 10, print_level = printlevel)
    rule_single_bidding = SDDP.DecisionRule(model_single_bidding; node = 1)
    sol_single_bidding = SDDP.evaluate(
        rule_single_bidding;
        incoming_state = Dict(Symbol("l[$(r.dischargepoint)]") => r.currentvolume for r in R),
        controls_to_record = [:x],
    )
    return 
end

function SingleOwnerScheduling(R::Array{Reservoir},
    K::Array{HydropowerPlant},
    y_initial::Dict{Int64, Float64},
    price::Vector{Float64},
    f_initial::Dict{Reservoir, Float64},
    Omega,
    P;
    mu_up = 1.0,
    mu_down = 0.1,
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    stopping_rule = [SDDP.BoundStalling(10, 1e-4)],
    T = 24,
    Stages = stages,
    optimizer = CPLEX.Optimizer)
    function subproblem_builder_single_short(subproblem::Model, node::Int64)
        # State Variables
        @variable(subproblem, 0 <= l[r = R] <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
        @variable(subproblem, ustart[k = K], Bin, SDDP.State, initial_value = 0)
        # Control Variables
        @variable(subproblem, y[t = 1:T])
        @variable(subproblem, Q[t = 1:T, r = R] >= 0)
        @variable(subproblem, 0 <= w[t = 1:T, k = K] <= k.spillreference * k.equivalent)
        @variable(subproblem, 0 <= u[t = 1:T, k = K], Bin)
        @variable(subproblem, 0 <= d[t = 1:T, k = K], Bin)
        @variable(subproblem, s[r = R] >= 0)
        @variable(subproblem, a[r = R])
        # Random Variables
        @variable(subproblem, f[r = R])
    
        # Transition function
        @constraint(subproblem, balance[r = R], l[r].out == l[r].in - sum(Q[t, r] for t in 1:T) + T * f[r] - s[r])
        @constraint(subproblem, planttrans1[k = K], ustart[k].in == u[1,k])
        @constraint(subproblem, planttrans2[k = K], ustart[k].out == u[T,k])
        # Constraints
        if node == 1
            @variable(subproblem, z_up[t = 1:T] >= 0)
            @variable(subproblem, z_down[t = 1:T] >= 0)
            @constraint(subproblem, obligation[t = 1:T], y[t] == sum(w[t,k] for k in K) + z_up[t] - z_down[t])
        end
        @constraint(subproblem, startup[t = 1:T-1, k = K], d[t,k] >= u[t+1,k] - u[t,k])
        @constraint(subproblem, activeplant[t = 1:T, k = K], w[t,k] <= u[t,k] * k.equivalent * k.spillreference)
        @constraint(subproblem, production[t = 1:T, k = K], w[t,k] <= sum(Q[t, r] for r in find_us_reservoir(k.reservoir)) * k.equivalent)
        @constraint(subproblem, watervalue[r = R], a[r] >= sum(k.equivalent for k in filter(k -> k.reservoir in find_ds_reservoirs(r), K)) * (l[r].in - l[r].out) * 0.5)
        # Parameterize Uncertainty
        SDDP.parameterize(subproblem, Omega, P) do om
            if node == 1
                for t in 1:T
                    JuMP.fix(y[t], y_initial[t])
                end
                for r in R
                    JuMP.fix(f[r], f_initial[r])
                end
                @stageobjective(subproblem, sum(price[t] * y[t] - mu_up * z_up[t] + mu_down * z_down[t]  - sum(S * d[t,k] for k in K) for t in 1:T) -sum(a[r] for r in R))
            else
                for r in R
                    JuMP.fix(f[r], om.inflow[r])
                end
                @stageobjective(subproblem, sum(om.price[t] * sum(w[t, k] for k in K) - sum(S * d[t, k] for k in K) for t in 1:T) - sum(a[r] for r in R))
            end
        end
    end
    model_single_short = SDDP.LinearPolicyGraph(subproblem_builder_single_short;
    stages = Stages, sense = :Max, upper_bound = Stages * T * sum(k.spillreference * k.equivalent for k in K), optimizer = optimizer)

    SDDP.train(model_single_short; iteration_limit = 10, stopping_rules = stopping_rule, print_level = printlevel)

    rule_single_short = SDDP.DecisionRule(model_single_short; node = 1)
    sol_single_short = SDDP.evaluate(rule_single_short; incoming_state = Dict(Symbol("l[$(r.dischargepoint)]") => r.currentvolume for r in R),
    controls_to_record = [:Q, :u, :d],)
    return
end

end

#=
 ------------------------------------------------------------------------ DEPRECATED --------------------------------------------------------------------------------------------------
 Sollte ich bald löschen, habe aber aktuell noch keinen Mut dazu.
# function Generation_Cuts(k::HydropowerPlant, n::Int64)
#     x = range(0,k.spillreference,n)
#     effs  = range(k.equivalent, k.equivalent/2, n)
#     y = x .* effs
#     # Get Slopes and y value of origin to obtain cut coefficients
#     f1 = [(y[i] - y[i-1])/(x[i] - x[i-1]) for i in 2:n] 
#     f2 = [y[i] - f1[i-1] * x[i] for i in 2:n]
#     return f1, f2, x, y
# end

# abstract type AbstractConfiguration end

# struct AnticipatoryConfig <: AbstractConfiguration end

# struct NonAnticipatoryConfig <: AbstractConfiguration end

# function add_state_variables(subproblem::Model, res::Vector{Reservoir}, res_real_initial::Dict{Reservoir, Float64},res_ind_initial::Dict{Reservoir, Float64}, ::AnticipatoryConfig)
#     @variables(subproblem, begin
#         0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
#         res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
#         Qnom[r = res], (SDDP.State, initial_value = 0)
#     end)
#     return
# end

# function add_state_variables(subproblem::Model, res::Vector{Reservoir}, res_real_initial::Dict{Reservoir, Float64}, res_ind_initial::Dict{Reservoir, Float64}, ::NonAnticipatoryConfig)
#     @variables(subproblem, begin
#         0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
#         res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
#         Qnom[r = res], (SDDP.State, initial_value = 0)
#     end)
#     return
# end

# function add_control_variables(subproblem::Model, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, plants_O::Vector{HydropowerPlant}, T::Int64, ::AnticipatoryConfig)
#     @variables(subproblem, begin
#         Qnom_change[r = res] >= 0
#         Qeff[k = plants_j, t = 1:T] >= 0
#         Qreal[r = res, t = 1:T] >= 0
#         Qadj[r = res] >= 0
#         P_Swap[r = res]
#         P_Over[k = plants_O] >= 0
#         BALANCE_INDICATOR[r = res], Bin
#     end)
#     return
# end

# function add_control_variables(subproblem::Model, res::Vector{Reservoir}, plants_j:: Vector{HydropowerPlant}, T::Int64, ::NonAnticipatoryConfig)
#     @variables(subproblem, begin
#         Qeff[k = plants_j, t = 1:T] >= 0
#         Qreal[r = res, t = 1:T] >= 0
#         Qnom_change[r = res] >= 0
#         BALANCE_INDICATOR[r = res], Bin
#     end)
#     return
# end

# function add_random_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, ::NonAnticipatoryConfig)
#     @variables(subproblem, begin
#         c[t = 1:T]
#         Qinflow[r = res] >= 0
#     end)
#     return
# end

# function add_random_variables(subproblem::Model, res::Vector{Reservoir}, T::Int64, ::AnticipatoryConfig)   
#     @variables(subproblem, begin
#         c[t = 1:T]
#         Qinflow[r = res] >= 0
#         Qnom_O[r = res]
#     end)
#     return
# end

# function add_stage_objective(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, j::Participant, mean_price::Dict{Reservoir, Float64}, T::Int64,::NonAnticipatoryConfig)
#     Qeff = subproblem[:Qeff]
#     c = subproblem[:c]
#     res_ind = subproblem[:res_ind]
#     if node == 1
#         @stageobjective(subproblem, 0)
#     else
#         @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
#         # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
#         + sum((res_ind[r].out - res_ind[r].in) * j.participationrate[r] * mean_price[r] for r in res))/1e3)
#     end
#     return
# end

# function add_stage_objective(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, j::Participant, mean_price::Dict{Reservoir, Float64}, T::Int64, ::AnticipatoryConfig)
#     Qeff = subproblem[:Qeff]
#     c = subproblem[:c]
#     res_ind = subproblem[:res_ind]
#     P_Swap = subproblem[:P_Swap]
#     if node == 1
#         # Objective function
#         @stageobjective(subproblem, 0)
#     else
#         # Objective Function
#         @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
#         + sum(c[t] * P_Swap[r] for t in 1:T for r in res)
#         # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
#         + sum((res_ind[r].out - res_ind[r].in)  * j.participationrate[r] * mean_price[r] for r in res))/1e3)
#     end
#     return
# end

# function add_transition_function(subproblem::Model, res::Vector{Reservoir}, node::Int64, T::Int64, Qref::Dict{Reservoir, Float64}, ::NonAnticipatoryConfig)
#     res_real = subproblem[:res_real]
#     res_ind = subproblem[:res_ind]
#     Qnom = subproblem[:Qnom]
#     Qnom_change = subproblem[:Qnom_change]
#     Qinflow = subproblem[:Qinflow]
#     if node == 1
#         for r in res
#             @constraint(subproblem, res_real[r].out == res_real[r].in)
#             @constraint(subproblem, res_ind[r].out == res_ind[r].in)
#             @constraint(subproblem, Qnom[r].out == Qnom_change[r])
#         end
#     else
#         for r in res
#             @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
#             @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
#             @constraint(subproblem, Qnom[r].out == Qnom_change[r])
#         end
#     end
#     return
# end

# function add_transition_function(subproblem::Model, res::Vector{Reservoir}, node::Int64, T::Int64, Qref::Dict{Reservoir, Float64}, ::AnticipatoryConfig)
#     res_real = subproblem[:res_real]
#     res_ind = subproblem[:res_ind]
#     Qnom = subproblem[:Qnom]
#     Qnom_change = subproblem[:Qnom_change]
#     Qinflow = subproblem[:Qinflow]
#     if node == 1
#         for r in res
#             @constraint(subproblem, res_real[r].out == res_real[r].in)
#             @constraint(subproblem, res_ind[r].out == res_ind[r].in)
#             @constraint(subproblem, Qnom[r].out == Qnom_change[r])
#         end
#     else
#         for r in res
#             @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
#             @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
#             @constraint(subproblem, Qnom[r].out == Qnom_change[r])
#         end
#     end
# end

# function add_stage_constraints(subproblem::Model, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, Qref::Dict{Reservoir, Float64}, BIG_M::Float64, T::Int64, stage_count::Int64, node::Int64, ::NonAnticipatoryConfig)
#     res_real = subproblem[:res_real]
#     res_ind = subproblem[:res_ind]
#     Qnom = subproblem[:Qnom]
#     BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
#     Qnom_change = subproblem[:Qnom_change]
#     Qeff = subproblem[:Qeff]
#     Qreal = subproblem[:Qreal]
#     if node == 1
#         for r in res
#             @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
#             @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
#             # Constraints
#             @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
#         end
#     else
#         for r in res
#             @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
#             @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qnom[r].in)
#             # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
#             # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
#             @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
#             @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
#         end
#         for k in plants_j
#             # @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
#             @constraint(subproblem, BALANCE_INDICATOR[k.reservoir] => {sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference})
#         end
#         for t in 1:T
#             for k in plants_j
#                 @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
#                 @constraint(subproblem, Qeff[k, t] <= k.spillreference)
#             end
#         end
#     end
#     return
# end

# function add_stage_constraints(subproblem::Model, node::Int64, res::Vector{Reservoir}, plants_j::Vector{HydropowerPlant}, plants_O::Vector{HydropowerPlant}, j::Participant, O::Participant, Qref::Dict{Reservoir, Float64}, BIG_M::Float64, T::Int64, stage_count::Int64, ::AnticipatoryConfig)
#     res_real = subproblem[:res_real]
#     res_ind = subproblem[:res_ind]
#     Qnom = subproblem[:Qnom]
#     BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
#     Qnom_change = subproblem[:Qnom_change]
#     Qeff = subproblem[:Qeff]
#     Qreal = subproblem[:Qreal]
#     Qadj = subproblem[:Qadj]
#     P_Over = subproblem[:P_Over]
#     P_Swap = subproblem[:P_Swap]
#     Qinflow = subproblem[:Qinflow]
#     Qnom_O = subproblem[:Qnom_O]
#     if node == 1
#         for k in plants_j
#             # @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
#             @constraint(subproblem, BALANCE_INDICATOR[k.reservoir] => {sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference})
#         end
#         for r in res
#             # Constraints
#             # @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
#             # @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
#             @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
#             @constraint(subproblem, BALANCE_INDICATOR[r] => {Qnom_change[r] <= Qref[r]})
#             @constraint(subproblem, !BALANCE_INDICATOR[r] => {0 <= res_ind[r].in})
#         end
#     else
#         for r in res
#             # Constraints
#             @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
#             @constraint(subproblem, Qadj[r] == (Qnom_O[r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
#             @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
#             @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
#             @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
#             @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
#         end
#         for k in plants_O
#             @constraint(subproblem, P_Over[k] >=  (sum(Qadj[r] for r in find_us_reservoir(k.reservoir)) - k.spillreference) * k.equivalent)
#             @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spillreference + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
#         end
#         for t in 1:T
#             for k in plants_j
#                 @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
#                 @constraint(subproblem, Qeff[k, t] <= k.spillreference)
#             end
#         end
#     end
#     return
# end

# function fix_random_variables(subproblem, res, node, T, Ω, P, ::NonAnticipatoryConfig)
#     c = subproblem[:c]
#     Qinflow = subproblem[:Qinflow]
#     if !(node == 1)
#         SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
#             for t in 1:T
#                 JuMP.fix(c[t], ω.price[t])
#             end
#             for r in res
#                 JuMP.fix(Qinflow[r], ω.inflow; force=true)
#             end
#         end
#     end
#     return
# end

# function fix_random_variables(subproblem, res, node, T, Ω, P, ::AnticipatoryConfig)
#     c = subproblem[:c]
#     Qinflow = subproblem[:Qinflow]
#     Qnom_O = subproblem[:Qnom_O]
#     if !(node == 1)
#         SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
#             for t in 1:T
#                 JuMP.fix(c[t], ω.price[t])
#             end
#             for r in res
#                 JuMP.fix(Qinflow[r], ω.inflow; force=true)
#                 try
#                     JuMP.fix(Qnom_O[r], ω.nomination[r], force=true)
#                 catch
#                     JuMP.fix(Qnom_O[r], ω.nomination[node][r], force=true)
#                 end
#             end
#         end 
#     end
#     return
# end

# function ShortTermOptimizationNoAnticipation(
#     all_res::Array{Reservoir},
#     j::Participant,
#     O::Participant,
#     plants_j::Array{HydropowerPlant},
#     iteration_count::Int64,
#     stage_count::Int64,
#     scenario_count::Int64,
#     T::Int64,
#     res_real_initial::Dict{Reservoir, Float64},
#     res_ind_initial::Dict{Reservoir, Float64},
#     Ω,
#     P,
#     Qref,
#     mean_price,
#     price_sample;
#     riskmeasure = SDDP.Expectation(),
#     printlevel = 1,
#     optimizer = CPLEX.Optimizer,
#     BIG_M = 5e4,
#     stall_bound = SDDP.BoundStalling(5, 1e-2),
#     config = NonAnticipatoryConfig()
#     )
#     res = filter(r -> j.participationrate[r] > 0, all_res)
#     function subproblem_builder(subproblem::Model, node::Int)
#         # Define State, Control and Random Variables
#         add_state_variables(subproblem, res, res_real_initial, res_ind_initial, config)
#         add_control_variables(subproblem, res, plants_j, T, config)
#         add_random_variables(subproblem, res, T, config)
#         # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
#         add_transition_function(subproblem, res, node, T, Qref, config)
#         # Add constraints for every stage
#         add_stage_constraints(subproblem, res, plants_j, Qref, BIG_M, T, stage_count, node, config)
#         # Fix Random Variables for nondeterministic stages
#         fix_random_variables(subproblem, res, node, T, Ω, P, config)
#         # Add Objective Function 
#         add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
#         return subproblem
#     end
#     # Define the policy graph structures and model
#     model = SDDP.LinearPolicyGraph(
#         subproblem_builder,
#         stages = stage_count,
#         sense = :Min,
#         lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
#         optimizer = optimizer
#     )
#     # Train the model
#     SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [stall_bound], print_level = printlevel, risk_measure = riskmeasure)
#     # obtain decision rule in all steps
#     rules = []
#     nominations = []
#     for node in 1:stage_count
#         rule = SDDP.DecisionRule(model; node = node)
#         solution = SDDP.evaluate(
#             rule;
#             incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
#             noise = (price = price_sample[rand(1:stage_count)][rand(1:scenario_count)], inflow = 0),
#             controls_to_record = [:Qnom]
#         )
#         Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
#         for r in filter(x -> !(x in res), all_res)
#             Qnom[r] = Qref[r]
#         end
#         append!(rules, [rule])
#         append!(nominations, [Qnom])
#     end
#     return model, rules, nominations
# end

# function ShortTermOptimizationAnticipation(
#     res::Array{Reservoir},
#     j::Participant,
#     O::Participant,
#     plants_O::Array{HydropowerPlant},
#     parts::Array{Participant},
#     plants_j::Array{HydropowerPlant},
#     ITERATION_COUNT::Int64,
#     stage_count::Int64,
#     SCENARIO_COUNT::Int64,
#     T::Int64,
#     res_real_initial::Dict{Reservoir, Float64},
#     res_ind_initial::Dict{Reservoir, Float64},
#     Ω,
#     P,
#     nom,
#     Qref,
#     mean_price,
#     price_sample;
#     riskmeasure = SDDP.Expectation(),
#     optimizer = CPLEX.Optimizer,
#     printlevel = 1,
#     BIG_M = 5e4,
#     stall_bound = SDDP.BoundStalling(5, 1e-2),
#     config = AnticipatoryConfig()
#     )
#     function subproblem_builder(subproblem::Model, node::Int)
#         # Add State, Control and Random Variables
#         add_state_variables(subproblem, res, res_real_initial, res_ind_initial, config)
#         add_control_variables(subproblem, res, plants_j, plants_O, T, config)
#         add_random_variables(subproblem, res, T, config)
#         # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
#         add_transition_function(subproblem, res, node, T, Qref, config)
#         # Add constraints for every stage
#         add_stage_constraints(subproblem, node, res, plants_j, plants_O, j, O, Qref, BIG_M, T, stage_count, config)
#         # Fix Random Variables for nondeterministic stages
#         fix_random_variables(subproblem, res, node, T, Ω, P, config)
#         # Add Objective Function 
#         add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
#         return subproblem
#     end
#     # define the policy graph structures and model
#     model = SDDP.LinearPolicyGraph(
#         subproblem_builder,
#         stages = stage_count,
#         sense = :Min,
#         lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
#         optimizer = optimizer
#     )
#     # train the model
#     SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound], risk_measure = riskmeasure, print_level = printlevel)
#     # obtain decision rule in all steps, as well as nominations for the reservoir situation.
#     rules = []
#     nominations = []
#     for node in 1:stage_count
#         rule = SDDP.DecisionRule(model; node = node)
#         solution = SDDP.evaluate(
#             rule;
#             incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
#             noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 0.1, nomination = nom[rand(1:SCENARIO_COUNT)]),
#             controls_to_record = [:Qnom]
#         )
#         Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
#         push!(rules, rule)
#         push!(nominations, Qnom)
#     end
#     return model, rules, nominations
# end

# function OptimizationAfterAdjustment(res::Array{Reservoir}, plants_j::Array{HydropowerPlant}, Qadj_par::Dict{Reservoir, Float64}, P_Swap_par::Dict{Reservoir, Float64}, c_par::Vector{Float64}, T::Int64; optimizer = CPLEX.Optimizer) 
#     model = Model(optimizer)
#     # Variables
#     @variable(model, Qreal[r = res, t=1:T] >= 0)
#     @variable(model, Qeff[k = plants_j,t=1:T] >= 0)
#     @variable(model, Qadj[r = res] >= 0)
#     @variable(model, P_Swap[r = res])
#     # Fix Variables which are given as parameters (because they are predetermined)
#     for r in res
#         JuMP.fix(Qadj[r], Qadj_par[r], force=true)
#         JuMP.fix(P_Swap[r], P_Swap_par[r])
#     end
#     for t in 1:T
#         JuMP.fix(c[t], c_par[t])
#     end
#     # Constraints -> Mostly production constraints determined thorugh spillage etc.
#     for r in res
#         @constraint(model, sum(Qreal[r, t]) == T * Qadj[r])
#     end
#     for t in 1:T
#         @constraint(model, sum(Qeff[k, t] for k in filter(x -> x.reservoir in find_ds_reservoirs(r), plants_j)) >= - P_Swap[r])
#         for k in plants_j
#             @constraint(model, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
#             @constraint(model, Qeff[k, t] <= k.spillreference)
#         end
#     end
#     # Objective -> Maximize revenue
#     @objective(model, Max, sum(c[t] * Qeff[k,t] for k in plants_j for t in 1:T) + sum(c[t] * P_Swap[r] for r in res for t in 1:T))
#     optimize!(model)
#     return Qeff, objective_value(model)
# end


# end

=#