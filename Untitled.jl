
vcf_filename = "test_4X_191.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
plot_types = "heatmap"
save_format = "html"
plot_title = "title_here"


using ViVa
ViVa.jupyter_main_new(vcf_filename,field_to_visualize,variant_filter,sample_filter,plot_types,save_format,plot_title,plot_labels)


using ViVa

