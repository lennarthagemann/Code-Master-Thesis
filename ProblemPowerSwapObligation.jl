using PlotlyJS, CSV, DataFrames, Distributions

df = DataFrame(hour = 1:24, PSwap = [10.0 for t in 1:24], Obligation = [25.0 + rand(Uniform(-5, 5)) for t in 1:24], Max = [35.0 for t in 1:24])
p = plot([
    bar(df, x=:hour, y=:PSwap, kind="bar", marker_color="yellow", name = "Power Swap"),
    bar(df, x=:hour, y=:Obligation, kind="bar", marker_color="lightblue", name = "Day Ahead Market Obligation"),
    bar(df, x=:hour, y=:Max, kind="scatter", mode="lines", marker_color="red",  name = "Maximum Production Capacity")
], Layout(barmode = "stack",
          xaxis_title = "Hour",
          yaxis_title = "Power"))

savefig(p, "C:\\Users\\lenna\\OneDrive - NTNU\\Master Thesis\\Final Presentation\\Images\\WarningExamplePowerSwap.png")