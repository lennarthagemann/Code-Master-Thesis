using JSON
using Plots
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation


filepath_Ljungan = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/TestDataWaterRegulation/Ljungan.json"
filepath_results = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Results/LambdaZero"
filename = "ResultsAnticipation4.json"
plot_savepath =" C://Users/lenna/OneDrive - NTNU/Master Thesis/Presentation VF/Images/ResultsAnticipation1.png"


res, plants, parts = read_data(filepath_Ljungan)
rounds = 7


json_dictionary_parts = Dict(
    parts[1] => "Name: Sydkraft\n",
    parts[2] => "Name: Fortum\n",
    parts[3] => "Name: Statkraft\n" )

json_dictionary_res = Dict(
    res[1] => "Reservoir with name: Flasjön\n",
    res[2] => "Reservoir with name: Holmsjön\n")

json_dictionary_plants = Dict(
    plants[1] => "A hydropowerplant with the name Flasjö and following information: \nName      : Flasjö\nReservoir : Flasjön\nEquivalent: 0.31\nSpill reference level: 0.58\n",
    plants[2] => "A hydropowerplant with the name Trangfors and following information: \nName      : Trangfors\nReservoir : Flasjön\nEquivalent: 0.72\nSpill reference level: 1.05\n",
    plants[3] => "A hydropowerplant with the name Rätan and following information: \nName      : Rätan\nReservoir : Flasjön\nEquivalent: 0.55\nSpill reference level: 1.15\n",
    plants[4] => "A hydropowerplant with the name Turinge and following information: \nName      : Turinge\nReservoir : Flasjön\nEquivalent: 0.19\nSpill reference level: 1.05\n",
    plants[5] => "A hydropowerplant with the name Bursnäs and following information: \nName      : Bursnäs\nReservoir : Flasjön\nEquivalent: 0.07\nSpill reference level: 1.05\n",
    plants[6] => "A hydropowerplant with the name Järnvägsforsen and following information: \nName      : Järnvägsforsen\nReservoir : Holmsjön\nEquivalent: 0.74\nSpill reference level: 1.45\n",
    plants[7] => "A hydropowerplant with the name Parteboda and following information: \nName      : Parteboda\nReservoir : Holmsjön\nEquivalent: 0.27\nSpill reference level: 1.4\n",
    plants[8] => "A hydropowerplant with the name Hermansboda and following information: \nName      : Hermansboda\nReservoir : Holmsjön\nEquivalent: 0.1\nSpill reference level: 1.4\n",
    plants[9] => "A hydropowerplant with the name Ljunga and following information: \nName      : Ljunga\nReservoir : Holmsjön\nEquivalent: 0.45\nSpill reference level: 1.45\n",
    plants[10] => "A hydropowerplant with the name Nederede and following information: \nName      : Nederede\nReservoir : Holmsjön\nEquivalent: 0.07\nSpill reference level: 2.4\n",
    plants[11] => "A hydropowerplant with the name Skallböle and following information: \nName      : Skallböle\nReservoir : Holmsjön\nEquivalent: 0.18\nSpill reference level: 2.65\n",
    plants[12] => "A hydropowerplant with the name Matfors and following information: \nName      : Matfors\nReservoir : Holmsjön\nEquivalent: 0.8\nSpill reference level: 2.5\n",
    plants[13] => "A hydropowerplant with the name Viforsen and following information: \nName      : Viforsen\nReservoir : Holmsjön\nEquivalent: 0.07\nSpill reference level: 1.6\n",)

function read_results(filepath_results, filename, json_dictionary_parts, json_dictionary_res, json_dictionary_plants, parts, res, plants, rounds)
    json_data = JSON.parsefile(filepath_results * "/" * filename)
    Qnoms =  [Dict((participant = p, reservoir = r) => json_data["Qnom"][i]["(participant = $(json_dictionary_parts[p]), reservoir = $(json_dictionary_res[r]))"] for p in parts for r in res) for i in 1:rounds]
    Qadjs = [Dict(r => json_data["Qadj"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    QadjTots = [Dict(r => json_data["QadjTot"][i][json_dictionary_res[r]] for r in res) for i in 1:rounds]
    P_Swaps = [Dict(p => Dict(r => json_data["P_Swap"][i][json_dictionary_parts[p]][json_dictionary_res[r]] for r in res) for p in parts) for i in 1:rounds]
    POvers = [Dict(p => Dict(plant => json_data["POver"][i][json_dictionary_parts[parts[1]]][json_dictionary_plants[plants[3]]] for plant in plants) for p in parts) for i in 1:rounds]
    ΣPOvers = [Dict(plant => json_data["ΣPOver"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    MaxEnergys = [Dict(plant => json_data["MaxEnergy"][i][json_dictionary_plants[plant]] for plant in plants) for i in 1:rounds]
    return Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys
end

Qnoms, Qadjs, QadjTots, P_Swaps, POvers, ΣPOvers, MaxEnergys = read_results(filepath_results, filename, json_dictionary_parts, json_dictionary_res, json_dictionary_plants, parts, res, plants, rounds)


function plotTotAdjustedFlow(rounds::Int64, parts::Array{Participant}, Qnoms, Qadjs; savepath = plot_savepath, save=true)
    x = 1:rounds
    y = []
    p = plot(x, [sum(Qadjs[i][r] for r in res) for i in 1:rounds], label = "Qadj", linewidth=3)
    plot!(x, [sum(Qnoms[i][(participant = parts[1], reservoir = r)] for r in res) for i in 1:rounds], label = "A")
    plot!(x, [sum(Qnoms[i][(participant = parts[2], reservoir = r)] for r in res) for i in 1:rounds], label = "B")
    plot!(x, [sum(Qnoms[i][(participant = parts[3], reservoir = r)] for r in res) for i in 1:rounds], label = "C")
    plot!(legend=:outerbottom, legendolumn=4)
    xlabel!("Auction Period")
    ylabel!("Qnom / Qadj")
    title!("(Total) Nominations and Adjusted Flow  ")
    if save == true
        savefig(savepath)
    end
    display(p)
end

    
plotTotAdjustedFlow(rounds, parts, Qnoms, Qadjs; save = true)
