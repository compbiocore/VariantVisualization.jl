
vcf_filename = "AC_gatk406_eh_PASS_withheader.vcf"
field_to_visualize = "genotype" 
variant_filter = "all"
sample_filter = "none" #"select_columns", "select_column_list.txt"
save_format = "pdf"
plot_title = "Jess_matrix"

vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only", "list", "significantList_for_proteinstructures.csv"#"range", "chr1:20000000-30000000"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
save_format = "pdf"
plot_title = "Example_2"




using ViVa

ViVa.jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)


using ViVa


vcf_tuple = ViVa.load_vcf(vcf_filename)
original_vcf = vcf_tuple[1]
df_vcf = vcf_tuple[2]
vcf = copy(original_vcf)
index = ViVa.format_reader(vcf, field_to_visualize)


function save_numerical_array1(x,y,z)

    df_withsamplenames = CSV.read(y, delim="\t", datarow = header_col+1, header = false, types=Dict(1=>String))
    samplenames=df_withsamplenames[1,10:size(df_withsamplenames,2)]
      samplenames=Matrix(samplenames)
      headings = hcat("chr","position")
      samplenames = hcat(headings,samplenames)
      chrlabels=z[:,1:2]

      chr_labeled_array_for_plotly=hcat(chrlabels, array_for_plotly)
      labeled_value_matrix_withsamplenames= vcat(samplenames,chr_labeled_array_for_plotly)

      writedlm("AC_gatk406_eh_PASS_withheader_value_matrix_.txt", labeled_value_matrix_withsamplenames, "\t")

end

value_matrix = ViVa.genotype_cell_searcher(vcf,index)
    array_for_plotly=value_matrix[:,10:size(value_matrix,2)]
    save_numerical_array1(array_for_plotly,vcf_filename,vcf)
