function f(x,y)
    return sqrt(x^3 + y^2)
end

# Easier: f(x,y) = sqrt(x^3 + y^2)
g(x) = log(x)
#=


for i in 1:10
    for j in 1:5
        println(f(i,j))
    end
end
=#
gf = g ∘ f # Komposition von Funktionen möglich

#println(2gf(10,5)) #Implizites Mulitiplizieren einer Zahl mit einer Variablen erlaubt.

# Ternäre Operatoren in Julia können Programme stark abkürzen
# <s1> && <s2> (short-circuit AND) Falls s1 falsch ist, wird false ausgegeben und s2 gar nicht mehr überprüft (AND kann nicht mehr wahr werden unabhängig von s2)
# <s1> || <s2> (short-circuit OR) Analog, hier wird falls s1 wahr ist s2 nicht mehr überprüft.
# x = (<s1> ? a : b ) Ternärer Operator, x wird auf a oder b gesetzt, je nachdem ob s1 wahr ist oder nicht.
c = 12
d = 10
b = true
x = b >= c ? "Das ist wahr." : "Das ist falsch."
println(x)
println(c > b && "Das ist wahr.")
println(d >= c && "Das ist wahr.")

# Matrixoperationen und elementweise
# p(x) = 4x + x^2
# A = [ [1,2,3][0,2,2][3,0,1] ]
# println(p(A)) #Matrixoperationen
# println(broadcast(p, A)) #Elementweise Operation
# println(p.(A)) #Analog, elementweise Operation

# Beliebige Anzahl an Funktionsparameter mittels ... (splat Operator)
collatz(el) = return (mod(el, 2) == 0  ? el/2 : 3el + 1)
function j(x...) 
    counter = []
    for el in x
        count = 0
        while el != 1
            el = collatz(el)
            count += 1
        end
        push!(counter, count)
    end
    return counter
end
# Mit ... als Suffix kann man nun beliebig viele
println(j(1,2,3,4,5))
println(j(12))

function symbolFn(n)
	subscriptn = String([Char(0x2080+d) for d in reverse(digits(n))])
	return Symbol("F",subscriptn)
end

println(symbolFn(123))

# Example of Expression and list comprehension

#Expression
begin
	F₁ = 1
	F₂ = 1
	for n ∈ 3:50
		eval(:( $(symbolFn(n)) = $(symbolFn(n-1)) + $(symbolFn(n-2))))
	end
end

# Print all 
[println( String(s), " = $(eval(:($s)))" ) for s ∈ symbolFn.(1:15)]

# Macros

macro multiple(n::Integer, expr) # with :: you can specify datatype, code will not run unless integer is given
    return quote
        $([esc(expr) for i ∈ 1:n]...)
    end
end

begin
	x = 1
	@multiple 10 x += 1
end;

#=
useful Macros:
@elapsed - Gibt die Zeit an, wie lange das Argument zur Bearbeitung braucht in Sekunden
@time - Wie @elapsed, speichert das Ergebnis zusätzlich in einer Variable
=#

begin
	inv([[1,0] [0,1]])
	@elapsed inv(rand(100,100))
end

#=
mit PrettyTables kann man Matrizen in schöne Tabellen umbauen.
=#

using PrettyTables

M = rand(4,4)
pretty_table(M)