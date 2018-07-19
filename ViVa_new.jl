
vcf_file_name = "variants.filtered.191_joint.vcf"
field_to_visualize = "read_depth" 
variant_filter = "all"
sample_filter = "select_columns", "select_column_list.txt"
save_format = "png"
plot_title = "Example_1"

vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only", "list", "significantList_for_proteinstructures.csv"#"range", "chr1:20000000-30000000"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
save_format = "pdf"
plot_title = "Example_2"




using ViVa

ViVa.jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)


readvcf = readlines(vcf_filename)

for row = 1:size(readvcf,1)

      if contains(readvcf[row], "#CHROM")
          header_col = row
          global header_col
          header_string = readvcf[row]
      end
end

df_withsamplenames = CSV.read(vcf_filename, delim="\t", datarow = header_col+1, header = false, types=Dict(1=>String))
Base.promote_rule(::Type{C}, ::Type{Any}) where {C <: CategoricalArrays.CatValue} = Any
matr=Matrix(df_withsamplenames)

#vcf_csv=CSV.read(vcf_filename, delim="\t", datarow = header_col+1, header = header_col, types=Dict(1=>String))

matr
