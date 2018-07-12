list = "DataFrames","CSV","PlotlyJS","Rsvg","Blink"
for i in list
Pkg.add(i)
end

Blink.AtomShell.install()
