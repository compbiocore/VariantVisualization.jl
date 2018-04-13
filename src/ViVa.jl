module ViVa

using DataFrames #use CSV.jl ? depwarnings
using PlotlyJS
using Rsvg
using Blink

const g_white = "400" #homo reference 0/0
const g_red = "800" #homo variant 1/1 1/2 2/2 1/3 2/3 3/3 4/4 5/5 6/6 etc
const g_pink = "600" #hetero variant 0/1 1/0 0/2 2/0 etc
const g_blue = "0" #no call ./.

include("vcf_utils.jl")
include("plot_utils.jl")

end # module
