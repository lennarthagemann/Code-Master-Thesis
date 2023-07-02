using JuMP, CPLEX

#Demonstriere ein Unit-Commitment Problem. Dies ist interessant für die Operation einer Wasserkraftanlage
# Wir möchten ein Produkt gewinnbringend verkaufen, dazu muss es kostengünstig produziert werden mit Generatoren.
# Die Generatoren können nur an oder aus sein, und es gilt über einen Zeitraum diese intelligent zu steuern.

c = [5.0, 4.0, 7.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 6.0] #Kostenvektor für Initalisierungskosten
g_min = [2.0, 4.0, 2.0, 0, 4.0, 0, 4.0, 2.0, 0, 0]
g_max = [10.0, 15.0, 20.0, 5.0, 11.0, 10.0, 15.0, 20.0, 5.0, 11.0]

N = length(c)
TIME_PERIODS = 5

model = Model(CPLEX.Optimizer)
set_silent(model)
@variable(model, x[1:N, 1:TIME_PERIODS], base_name = "Unit Status", Bin)  # Binary variable for status of generator
@objective(model, Min, sum(c' * x[t])) # Minimize costs over all time steps
