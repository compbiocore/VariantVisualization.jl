#this is what user enters in first cell before running whole notebook.
vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only", "list", "significantList_for_proteinstructures.csv"#"range", "chr1:20000000-30000000"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
# select_columns must have "sample_id_list.t", "reorder_columns", must have "phenotype_matrix.csv" & "phenotype_label_to_sort_by
save_format = "png"
plot_title = "title_here"



using DataFrames #use CSV.jl ? depwarnings
using CSV
using PlotlyJS
using Rsvg
using Blink
using ViVa

function jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)

#A)Load vcf and identify field index
vcf_tuple = ViVa.load_vcf(vcf_filename)
original_vcf = vcf_tuple[1]
df_vcf = vcf_tuple[2]

index = ViVa.format_reader(original_vcf, field_to_visualize)

#retain original version of vcf for second run
vcf = copy(original_vcf)

#B)Apply filters # issue - if vector need to use in, if string need to use contains() / sample_filters may be both - how to conditionally use correct function

#1)Variant filters
for i = 1:size(variant_filter,1)
    if variant_filter[i] == "pass_only"
        println("selecting pass_only variants")
        vcf=vcf[(vcf[:,7].== "PASS"),:]

    elseif variant_filter[i] == "range"
        chr_range = variant_filter[i+1]
        println("selecting variants within $chr_range")
        vcf = ViVa.chromosome_range_vcf_filter(chr_range,vcf)

    elseif variant_filter[i] == "list"
        siglist_file = variant_filter[i+1]

        siglist = (load_siglist(siglist_file))

        vcf = ViVa.sig_list_vcf_filter(vcf,siglist)

    end
end

#2) Sample Filters
for i = 1:size(sample_filter,1)
    if sample_filter[i] == "reorder_columns"
        list =sample_filter[i+1]
        key = sample_filter[i+2]
        println("sorting columns by $key in $list")
        vcf = load_sort_phenotype_matrix(list, key, vcf, df_vcf)

    elseif sample_filter[i] == "select_columns"
        id_list = sample_filter[i+1]
        println("selecting columns to match $id_list")
        vcf = ViVa.select_columns(id_list, vcf, df_vcf)
    end
end

#C) Convert filtered vcf to value matrix then plot and save
if field_to_visualize == "genotype"

    value_matrix = ViVa.genotype_cell_searcher(vcf,index)
    array_for_plotly=value_matrix[:,10:size(value_matrix,2)]
    graphic = ViVa.genotype_heatmap2(array_for_plotly,plot_title)
    PlotlyJS.savefig(graphic, "$plot_title.$save_format")

elseif field_to_visualize == "read_depth"

    value_matrix=ViVa.dp_cell_searcher(vcf,index)
    array_for_plotly=value_matrix[:,10:size(value_matrix,2)]
    graphic = ViVa.dp_heatmap2(array_for_plotly,plot_title)
    PlotlyJS.savefig(graphic, "$plot_title.$save_format")

end

end #end of jupyter_main()

jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)
