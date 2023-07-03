module WaterRegulation

using JSON
using Printf
using JuMP
using CPLEX
using SDDP
export HydropowerPlant, Reservoir, Participant, adjust_flow!, calculate_balance, update_reservoir, update_ind_reservoir, Calculate_Ersmax, Calculate_POver, power_swap, find_us_reservoir, find_ds_reservoirs, connect_reservoirs, read_nomination, read_data, water_regulation, OtherParticipant, CalculateQmax, Calculate_Qover, partAvg, SimplePartAvg, SumPartAvg, calculate_produced_power, total_power, ShortTermOptimizationNoAnticipation, ShortTermOptimizationAnticipation
 

mutable struct Reservoir
    dischargepoint::String
    totalvolume::Float64
    currentvolume::Float64
    maxvolume::Float64
    upstream_reservoir::Union{Array{Reservoir}, Nothing}
    downstream_reservoir::Union{Reservoir, Nothing}
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
        if typeof(dr.upstream_reservoir) == Nothing
            dr.upstream_reservoir = [r]
        else
            push!(dr.upstream_reservoir, r)
        end
        r.downstream_reservoir = dr
    end
end

"""
Find all upstream reservoirs of a given discharge point, including itself. Relevant for the calculation of any participation rate,
as participants can nominate for discharge at any upstream reservoir of their powerplant.
"""
function find_us_reservoir(r::Reservoir)::Array{Reservoir}
    queue::Array{Reservoir} = [r]
    us_res::Array{Reservoir} = []
    # Write while loop that terminates when queue is empty. Add upstream_reservoirs of elements in queue during iteration.
    while !isempty(queue)
        current_res::Reservoir = popfirst!(queue)
        @assert (typeof(current_res) == Reservoir)
        push!(us_res, current_res)
        if typeof(current_res.upstream_reservoir) != Nothing
            for res in current_res.upstream_reservoir
                push!(queue, res)
            end 
        end
    end
    return us_res
end

function find_ds_reservoirs(r::Reservoir)::Array{Reservoir}
    queue::Array{Reservoir} = [r]
    ds_res::Array{Reservoir} = []
    while !isempty(queue)
        current_res::Reservoir = popfirst!(queue)
        @assert (typeof(current_res) == Reservoir)
        push!(ds_res, current_res)
        if typeof(current_res.downstream_reservoir) != Nothing
            push!(queue, current_res.downstream_reservoir)
        end
    end
    return ds_res
end

struct HydropowerPlant
    name::String
    reservoir::Reservoir
    equivalent::Float64
    spill_reference_level::Float64
    function HydropowerPlant(name::String, reservoir::Reservoir, equivalent::Float64, spill_reference_level::Float64)
        return new(name, reservoir, equivalent, spill_reference_level)
    end
    function HydropowerPlant()
        return new("TestPlant", Reservoir(), 0.1, 1)
    end
end

struct Participant
    name::String
    plants::Array{HydropowerPlant}
    participationrate::Dict{Reservoir, Float64} 
    individual_reservoir::Dict{Reservoir, Float64}
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
    println(io,"Power Plant: ",  hp.name)
end

function Base.show(io::IO, r::Reservoir)
    println(io, "Reservoir with name: ", r.dischargepoint)
end

function Base.show(io::IO, p::Participant)
    println(io, "Name: ", p.name)
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
        nom.part.individual_reservoir[nom.res] += (Qref[nom.res] - value)*24
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
#    Qnom_cleaned = adjusted_nominations(Qnom)
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
                if Qadj[r] > Qmaxr[r] && !(nom.participant.individual_reservoir[nom.reservoir]< 0) 
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
function update_reservoir(Qadj_All::Dict{Reservoir, Float64}, T::Int64)
   for (d, Qadj) in Qadj_All
        if haskey(Qadj_All, d.upstream_reservoir)
            d.currentvolume -= Qadj*T + sum([adj for adj in Qadj_All[d.upstream_reservoir]])
        else
            d.currentvolume -= Qadj*T
        end
    end
end

"""
Update the individual reservoir by updating the balance of every participant through their new nomination.
Every participant has a field balance, it is the updated by the difference of nomination value and reference flow at every reservoir.
(Times 24, as it is done for the entire day)
"""
function update_ind_reservoir(Qnom::Dict{NamedTuple{(:participant, :reservoir), Tuple{Participant, Reservoir}}, Float64}, Qref::Dict{Reservoir, Float64})
    for (nom, value) in Qnom
        if haskey(nom.participant.individual_reservoir, nom.reservoir)
            nom.participant.individual_reservoir[nom.reservoir] += (Qref[nom.reservoir] - value)*24
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
    plants = [HydropowerPlant(p["name"], resdict[p["reservoir"]], float(p["equivalent"]), float(p["spill_reference_level"])) for p in d["HydropowerPlants"]]
    # Unpack Participants
    plantdict = Dict(p.name => p for p in plants)
    parts = [Participant(part["name"], [plantdict[pname] for pname in part["plants"]],  Dict(resdict[k] => float(v) for (k,v) in part["individual_reservoir"]), res) for part in d["Participants"]]
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
                Qmax[p] =  max(min(p.spill_reference_level, Qref[p.reservoir]), value)
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
        Ersmax[k] = max((QadjTot[k.reservoir] - k.spill_reference_level) * k.equivalent, 0)
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
Generalize code by adding configurations:
AnticipatoryConfiguration -> add variables, constraints, etc. relevant to anticipatory problem specification
NonAnticipatoryConfiguration -> likewise.
"""
abstract type AbstractConfiguration end

struct AnticipatoryConfig <: AbstractConfiguration end

struct NonAnticipatoryConfig <: AbstractConfiguration end

function add_state_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    ::AnticipatoryConfig
    )
    @variables(subproblem, begin
        0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
        res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
        Qnom[r = res], (SDDP.State, initial_value = 0)
    end)
    return
end

function add_state_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    ::NonAnticipatoryConfig
    )
    @variables(subproblem, begin
        0 <= res_real[r = res] <= r.maxvolume, (SDDP.State, initial_value = res_real_initial[r])
        res_ind[r = res], (SDDP.State, initial_value = res_ind_initial[r])
        Qnom[r = res], (SDDP.State, initial_value = 0)
    end)
    return
end

function add_control_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    plants_j::Vector{HydropowerPlant},
    plants_O::Vector{HydropowerPlant},
    T::Int64,
    ::AnticipatoryConfig
    )
    @variables(subproblem, begin
        Qnom_change[r = res] >= 0
        Qeff[k = plants_j, t = 1:T] >= 0
        Qreal[r = res, t = 1:T] >= 0
        Qadj[r = res] >= 0
        P_Swap[r = res]
        P_Over[k = plants_O] >= 0
        BALANCE_INDICATOR[r = res], Bin
    end)
    return
end

function add_control_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    plants_j:: Vector{HydropowerPlant},
    T::Int64,
    ::NonAnticipatoryConfig
    )
    @variables(subproblem, begin
        Qeff[k = plants_j, t = 1:T] >= 0
        Qreal[r = res, t = 1:T] >= 0
        Qnom_change[r = res] >= 0
        BALANCE_INDICATOR[r = res], Bin
    end)
    return
end

function add_random_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    T::Int64,
    ::NonAnticipatoryConfig
    )
    @variables(subproblem, begin
        c[t = 1:T]
        Qinflow[r = res] >= 0
    end)
    return
end

function add_random_variables(
    subproblem::Model,
    res::Vector{Reservoir},
    T::Int64,
    ::AnticipatoryConfig
    )   
    @variables(subproblem, begin
        c[t = 1:T]
        Qinflow[r = res] >= 0
        Qnom_O[r = res]
    end)
    return
end

function add_stage_objective(
    subproblem::Model,
    node::Int64,
    res::Vector{Reservoir},
    plants_j::Vector{HydropowerPlant},
    j::Participant,
    mean_price::Dict{Reservoir, Float64},
    T::Int64,
    ::NonAnticipatoryConfig
    )
    Qeff = subproblem[:Qeff]
    c = subproblem[:c]
    res_ind = subproblem[:res_ind]
    if node == 1
        @stageobjective(subproblem, 0)
    else
        @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
        # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
        + sum((res_ind[r].out - res_ind[r].in) * j.participationrate[r] * mean_price[r] for r in res))/1e3)
    end
    return
end

function add_stage_objective(
    subproblem::Model,
    node::Int64,
    res::Vector{Reservoir},
    plants_j::Vector{HydropowerPlant},
    j::Participant,
    mean_price::Dict{Reservoir, Float64},
    T::Int64,
    ::AnticipatoryConfig
    )
    Qeff = subproblem[:Qeff]
    c = subproblem[:c]
    res_ind = subproblem[:res_ind]
    P_Swap = subproblem[:P_Swap]
    if node == 1
        # Objective function
        @stageobjective(subproblem, 0)
    else
        # Objective Function
        @stageobjective(subproblem,  -(sum(c[t] * Qeff[k, t] * k.equivalent for t in 1:T for k in plants_j)
        + sum(c[t] * P_Swap[r] for t in 1:T for r in res)
        # + sum(((j.participationrate[r])/(j.participationrate[r] + O.participationrate[r])) * (res_real[r].out - res_real[r].in)  * j.participationrate[r] * mean_price[r] for r in res) 
        + sum((res_ind[r].out - res_ind[r].in)  * j.participationrate[r] * mean_price[r] for r in res))/1e3)
    end
    return
end

function add_transition_function(
    subproblem::Model,
    res::Vector{Reservoir},
    node::Int64,
    T::Int64,
    Qref::Dict{Reservoir, Float64},
    ::NonAnticipatoryConfig
    )
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    Qnom_change = subproblem[:Qnom_change]
    Qinflow = subproblem[:Qinflow]
    if node == 1
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in)
            @constraint(subproblem, res_ind[r].out == res_ind[r].in)
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
        end
    else
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
            @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
        end
    end
    return
end

function add_transition_function(
    subproblem::Model,
    res::Vector{Reservoir},
    node::Int64,
    T::Int64,
    Qref::Dict{Reservoir, Float64},
    ::AnticipatoryConfig
    )
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    Qnom_change = subproblem[:Qnom_change]
    Qinflow = subproblem[:Qinflow]
    if node == 1
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in)
            @constraint(subproblem, res_ind[r].out == res_ind[r].in)
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
        end
    else
        for r in res
            @constraint(subproblem, res_real[r].out == res_real[r].in - T * (Qnom[r].in - Qinflow[r]))
            @constraint(subproblem, res_ind[r].out == res_ind[r].in - T * (Qnom[r].in - Qref[r]))
            @constraint(subproblem, Qnom[r].out == Qnom_change[r])
        end
    end
end

function add_stage_constraints(
    subproblem::Model,
    res::Vector{Reservoir},
    plants_j::Vector{HydropowerPlant},
    Qref::Dict{Reservoir, Float64},
    BIG_M::Float64,
    T::Int64,
    stage_count::Int64,
    node::Int64,
    ::NonAnticipatoryConfig
    )
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
    Qnom_change = subproblem[:Qnom_change]
    Qeff = subproblem[:Qeff]
    Qreal = subproblem[:Qreal]
    if node == 1
        for r in res
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            for k in plants_j
                @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
            end
            # Constraints
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
        end
    else
        for r in res
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qnom[r].in)
            @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
        end
        for k in plants_j
            @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
        end
        for t in 1:T
            for k in plants_j
                @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
        end
    end
    return
end

function add_stage_constraints(
    subproblem::Model,
    node::Int64,
    res::Vector{Reservoir},
    plants_j::Vector{HydropowerPlant},
    plants_O::Vector{HydropowerPlant},
    j::Participant,
    O::Participant,
    Qref::Dict{Reservoir, Float64},
    BIG_M::Float64,
    T::Int64,
    stage_count::Int64,
    ::AnticipatoryConfig
    )
    res_real = subproblem[:res_real]
    res_ind = subproblem[:res_ind]
    Qnom = subproblem[:Qnom]
    BALANCE_INDICATOR = subproblem[:BALANCE_INDICATOR]
    Qnom_change = subproblem[:Qnom_change]
    Qeff = subproblem[:Qeff]
    Qreal = subproblem[:Qreal]
    Qadj = subproblem[:Qadj]
    P_Over = subproblem[:P_Over]
    P_Swap = subproblem[:P_Swap]
    Qinflow = subproblem[:Qinflow]
    Qnom_O = subproblem[:Qnom_O]
    if node == 1
        for k in plants_j
            @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
        end
        for r in res
            # Constraints
            @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
        end
    else
        for r in res
            # Constraints
            @constraint(subproblem, stage_count * T * Qnom_change[r] <= res_real[r].out)
            @constraint(subproblem, Qadj[r] == (Qnom_O[r] * O.participationrate[r] + Qnom[r].in * j.participationrate[r])/(O.participationrate[r] + j.participationrate[r]))
            @constraint(subproblem, P_Swap[r] ==  (Qnom[r].in - Qadj[r]) * j.participationrate[r] - sum(P_Over[k] for k in filter(k -> k.reservoir == r, plants_O)))
            @constraint(subproblem, sum(Qreal[r, t] for t in 1:T) == T * Qadj[r])
            @constraint(subproblem, Qnom_change[r] <= Qref[r] + BIG_M *(1 - BALANCE_INDICATOR[r]))
            @constraint(subproblem, 0 <= res_ind[r].in + BIG_M * BALANCE_INDICATOR[r])
        end
        for k in plants_O
            @constraint(subproblem, P_Over[k] >=  (sum(Qadj[r] for r in find_us_reservoir(k.reservoir)) - k.spill_reference_level) * k.equivalent)
            @constraint(subproblem, sum(Qnom_change[r_up] for r_up in find_us_reservoir(k.reservoir)) <= k.spill_reference_level + BIG_M * (1 - BALANCE_INDICATOR[k.reservoir]))
        end
        for t in 1:T
            for k in plants_j
                @constraint(subproblem, Qeff[k, t] <= sum(Qreal[r_up, t] for r_up in find_us_reservoir(k.reservoir)))
                @constraint(subproblem, Qeff[k, t] <= k.spill_reference_level)
            end
        end
    end
    return
end

function fix_random_variables(subproblem, res, node, T, Ω, P, ::NonAnticipatoryConfig)
    c = subproblem[:c]
    Qinflow = subproblem[:Qinflow]
    if !(node == 1)
        SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
            for t in 1:T
                JuMP.fix(c[t], ω.price[t])
            end
            for r in res
                JuMP.fix(Qinflow[r], ω.inflow; force=true)
            end
        end
    end
    return
end

function fix_random_variables(subproblem, res, node, T, Ω, P, ::AnticipatoryConfig)
    c = subproblem[:c]
    Qinflow = subproblem[:Qinflow]
    Qnom_O = subproblem[:Qnom_O]
    if !(node == 1)
        SDDP.parameterize(subproblem, Ω[node], P[node]) do ω
            for t in 1:T
                JuMP.fix(c[t], ω.price[t])
            end
            for r in res
                JuMP.fix(Qinflow[r], ω.inflow; force=true)
                try
                    JuMP.fix(Qnom_O[r], ω.nomination[r], force=true)
                catch
                    JuMP.fix(Qnom_O[r], ω.nomination[node][r], force=true)
                end
            end
        end 
    end
    return
end

function ShortTermOptimizationNoAnticipation(
    all_res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_j::Array{HydropowerPlant},
    iteration_count::Int64,
    stage_count::Int64,
    scenario_count::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Ω,
    P,
    Qref,
    mean_price,
    price_sample;
    riskmeasure = SDDP.Expectation(),
    printlevel = 1,
    optimizer = CPLEX.Optimizer,
    BIG_M = 1e5,
    stall_bound = SDDP.BoundStalling(5, 1e-2),
    config = NonAnticipatoryConfig()
    )
    res = filter(r -> j.participationrate[r] > 0, all_res)
    function subproblem_builder(subproblem::Model, node::Int)
        # Define State, Control and Random Variables
        add_state_variables(subproblem, res, res_real_initial, res_ind_initial, config)
        add_control_variables(subproblem, res, plants_j, T, config)
        add_random_variables(subproblem, res, T, config)
        # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
        add_transition_function(subproblem, res, node, T, Qref, config)
        # Add constraints for every stage
        add_stage_constraints(subproblem, res, plants_j, Qref, BIG_M, T, stage_count, node, config)
        # Fix Random Variables for nondeterministic stages
        fix_random_variables(subproblem, res, node, T, Ω, P, config)
        # Add Objective Function 
        add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
        return subproblem
    end
    # Define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
        optimizer = optimizer
    )
    # Train the model
    SDDP.train(model; iteration_limit = iteration_count, stopping_rules = [stall_bound], print_level = printlevel, risk_measure = riskmeasure)
    # obtain decision rule in all steps
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:scenario_count)], inflow = 0),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        for r in filter(x -> !(x in res), all_res)
            Qnom[r] = Qref[r]
        end
        append!(rules, [rule])
        append!(nominations, [Qnom])
    end
    return model, rules, nominations
end


function ShortTermOptimizationAnticipation(
    res::Array{Reservoir},
    j::Participant,
    O::Participant,
    plants_O::Array{HydropowerPlant},
    parts::Array{Participant},
    plants_j::Array{HydropowerPlant},
    ITERATION_COUNT::Int64,
    stage_count::Int64,
    SCENARIO_COUNT::Int64,
    T::Int64,
    res_real_initial::Dict{Reservoir, Float64},
    res_ind_initial::Dict{Reservoir, Float64},
    Ω,
    P,
    nom,
    Qref,
    mean_price,
    price_sample;
    riskmeasure = SDDP.Expectation(),
    optimizer = CPLEX.Optimizer,
    printlevel = 1,
    BIG_M = 1e5,
    stall_bound = SDDP.BoundStalling(5, 1e-2),
    config = AnticipatoryConfig()
    )
    function subproblem_builder(subproblem::Model, node::Int)
        # Add State, Control and Random Variables
        add_state_variables(subproblem, res, res_real_initial, res_ind_initial, config)
        add_control_variables(subproblem, res, plants_j, plants_O, T, config)
        add_random_variables(subproblem, res, T, config)
        # Add Transition Function for reservoir levels, individual reservoirs and propagation of nomination to next stage
        add_transition_function(subproblem, res, node, T, Qref, config)
        # Add constraints for every stage
        add_stage_constraints(subproblem, node, res, plants_j, plants_O, j, O, Qref, BIG_M, T, stage_count, config)
        # Fix Random Variables for nondeterministic stages
        fix_random_variables(subproblem, res, node, T, Ω, P, config)
        # Add Objective Function 
        add_stage_objective(subproblem, node, res, plants_j, j, mean_price, T, config)
        return subproblem
    end
    # define the policy graph structures and model
    model = SDDP.LinearPolicyGraph(
        subproblem_builder,
        stages = stage_count,
        sense = :Min,
        lower_bound = -sum(r.maxvolume * mean_price[r] * j.participationrate[r] for r in res)/1e3,
        optimizer = optimizer
    )
    # train the model
    SDDP.train(model; iteration_limit = ITERATION_COUNT, stopping_rules = [stall_bound], risk_measure = riskmeasure, print_level = printlevel)
    # obtain decision rule in all steps, as well as nominations for the reservoir situation.
    rules = []
    nominations = []
    for node in 1:stage_count
        rule = SDDP.DecisionRule(model; node = node)
        solution = SDDP.evaluate(
            rule;
            incoming_state = merge(Dict(Symbol("res_real[$(r)]") => res_real_initial[r] for r in res), Dict(Symbol("res_ind[$(r)]") => res_ind_initial[r] for r in res)),
            noise = (price = price_sample[rand(1:stage_count)][rand(1:SCENARIO_COUNT)], inflow = 0.1, nomination = nom[rand(1:SCENARIO_COUNT)]),
            controls_to_record = [:Qnom]
        )
        Qnom = Dict(r => solution.controls[:Qnom][r].out for r in res)
        push!(rules, rule)
        push!(nominations, Qnom)
    end
    return model, rules, nominations
end

end