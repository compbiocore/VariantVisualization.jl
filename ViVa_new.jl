
vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only", "list", "significantList_for_proteinstructures.csv"#"range", "chr1:20000000-30000000"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
save_format = "pdf"
plot_title = "title_here"
;




using DataFrames #use CSV.jl ? depwarnings
using CSV
using PlotlyJS
using Rsvg
using Blink
using ViVa

@time ViVa.jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)
