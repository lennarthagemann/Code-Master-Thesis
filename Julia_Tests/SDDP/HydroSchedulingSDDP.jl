using JuMP
using CPLEX
using Distributions
using LinearAlgebra
using SDDP
try
    using Revise
catch e
    @warn "Error initializing Revise" exception=(e, catch_backtrace())
end
includet("C:/Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/WaterRegulation.jl")
using .WaterRegulation

filepath_systemA = "C://Users/lenna/OneDrive - NTNU/Code Master Thesis/Water_Regulation/SimpleReservoirSystem.json"
res, plants, parts = read_data(filepath_systemA)

r = res[1]
p = parts[1]
Pmax = 10000
NODE_COUNT = 52
pplants = filter(x -> x.reservoir == r, p.plants)

graph = SDDP.UnicyclicGraph(0.8; num_nodes = NODE_COUNT)

function subproblem_builder(subproblem::Model, node::Int)
    @variable(subproblem, 0 <= volume <= r.maxvolume, SDDP.State, initial_value = r.currentvolume)
    @variables(subproblem, begin
        0 <= Q
        0 <= P[1:length(pplants)] 
    end)
    # Random Variables
    @variable(subproblem, inflow)
    Ω = [0.0, 500.0, 1000.0]
    P = [1 / 3, 1 / 3, 1 / 3]
    SDDP.parameterize(subproblem, Ω, P) do ω
        return JuMP.fix(inflow, ω)
    end
    @constraints(
        subproblem,
        begin
            volume.out == volume.in - Q + inflow
        end
    )
    for k in eachindex(pplants)
        @constraints(
            subproblem,
            begin
                P[k] <= Q * k.equivalent
                P[k] <= pplants[k].spill_reference_level * pplants[k].equivalent
            end
        )
    end
    prices = rand(Uniform(0, 100), NODE_COUNT)
    @stageobjective(subproblem, prices[node] * P)
    return subproblem
end

model = SDDP.PolicyGraph(
    subproblem_builder,
    graph;
    sense = :Max,
    upper_bound = 100 * NODE_COUNT * r.maxvolume ,
    lower_bound = 0,
    optimizer = CPLEX.Optimizer
)

SDDP.train(model; iteration_limit = 10)