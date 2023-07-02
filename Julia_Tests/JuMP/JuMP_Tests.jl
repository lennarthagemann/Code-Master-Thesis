using JuMP, HiGHS, Plots, PrettyTables

#1. Wichtige Features aus Julia (noch nicht in andere Tests abgedeckt)

x = 25
# println("Der Wert von x ist $(x)") # Mit $() kann man f-Strings aus Python nutzen
# println("Der Wert von x ist $(eval(:x))") # Mit : können wir das Symbol der Variable genau ansprechen, Datenstruktur des Compilers das Symbol eindutig identifiziert. Symbole brauchen weniger speicher

#Man kann eigene Datenstrukturen definieren, auch wenn Julia keine objektorientierte Sprache per se ist.
# Die Attribute sind mutable, es sei denn wir definieren es als 'mutable struct'.
# !Man sollte Datenstrukturen in eigene Module packen, da diese sonst sich nicht mehr verändern lassen!
module HydropowerPlants
export HydropowerPlant
struct HydropowerPlant
    ResSize::Int
    equiv::Float32
    Turbines::Int
    Name::String
    River::String
end
end

hp = HydropowerPlants.HydropowerPlant(12340000, 0.7, 4, "Vattenfall", "Glomma")

#Anonyme Funktionen: Funktioniere ähnlich wie in Java, nur etwas leichter und weniger boilerplate

a = [i for i in 1:5]
# println(map(x -> x^2, a))

#2. Einführung in JuMP

model = Model(HiGHS.Optimizer) #Gebe einem Optimierungsmodell einen Optimierer, der dieses nachher lösen soll. Model Objekte werden sukzessive spezifiziert.
# set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-8) #Setzt den Abbruchparameter der Opimierung auf 10^-8. Einstellungen spezifisch zum Solver, hier CPLEX

@variable(model, 0 <= x <= 2.5, Int)
@variable(model, 0 <= y <= 2.5, Int)
@objective(model, Max, y) #Zielfunktion die minimiert/maximiert werden soll

@constraint(model, c1, y - x <= 1) #Constraints werden einem model zugeordnet und bekommen einen Namen, enthalten schließlich die mathematische Restrikition als Expression

optimize!(model) # Nach Aufruf wird das Problem gelöst (! am Ende hat keine besondere Bedeutung, Konvention das Funktionen die die Argumente mutieren können mit ! aufhören.)

println("x : $(value.(x))")
println("y : $(value.(y))")

#Wir können mehrere Variablen mittels Matrizen zu einem Model hinzufügen. Diese haben drei Typen: Arrays, DenseAxisArrays, SparseAxisArrays

model2  = Model(HiGHS.Optimizer)
# set_optimizer_attribute(model2, "CPX_PARAM_EPINT", 1e-16)

@variable(model2, a[1:5,1:5])

@objective(model2, Max, sum(a))

@constraint(model2,[i in 1:5], a[i,:] .<= 1)

optimize!(model2)