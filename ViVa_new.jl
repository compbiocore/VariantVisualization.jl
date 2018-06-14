
vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "read_depth" 
variant_filter = "all"
sample_filter = "select_columns", "select_column_list.txt"
save_format = "png"
plot_title = "Example_1"
;

using DataFrames #use CSV.jl ? depwarnings
using CSV
using PlotlyJS
using Rsvg
using Blink
using ViVa

ViVa.jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)

