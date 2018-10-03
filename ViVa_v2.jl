
vcf_filename = "test_4X_191.vcf"
field_to_visualize = "read_depth" # field(s) = "read_depth", *etc.
variant_filter = ["range","chr4:0-4000000000"] #"range","chr4:0-4000000000"
sample_filter = ["reorder_columns", "sample_phenotype_matrix.csv", "case_control_status"] # "select_columns", "select_column_list.txt",
plot_types = "sample_line_chart"
save_format = "html"
plot_title = "Read_depth_test"


using ViVa
using GeneticVariation
using VCFTools
plot=ViVa.jupyter_main_new(vcf_filename,field_to_visualize,variant_filter,sample_filter,plot_types,save_format,plot_title)

