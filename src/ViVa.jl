module ViVa

using DataFrames #use CSV.jl ? depwarnings
using PlotlyJS
using Rsvg
using Blink
using CSV
using VCFTools

#include("vcf_utils.jl")
#include("plot_utils.jl")

#end #module

#=
const g_white = "400" #homo reference 0/0
const g_red = "800" #homo variant 1/1 1/2 2/2 1/3 2/3 3/3 4/4 5/5 6/6 etc
const g_pink = "600" #hetero variant 0/1 1/0 0/2 2/0 etc
const g_blue = "0" #no call ./.
=#

export
    format_reader,
    load_vcf,
    clean_column1!,
    genotype_cell_searcher_maf_correction,
    genotype_cell_searcher,
    dp_cell_searcher,
    load_siglist,
    sig_list_vcf_filter,
    chromosome_range_vcf_filter,
    load_sort_phenotype_matrix,
    reorder_columns,
    select_columns,
    genotype_heatmap2,
    dp_heatmap2,
    jupyter_main,
    save_numerical_array

#include("vcf_utils.jl")
include("type_defined_vcf_utils.jl")
include("plot_utils.jl")
include("notebook_utils.jl")

end # module
