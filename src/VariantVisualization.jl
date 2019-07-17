module VariantVisualization

using DataFrames
using PlotlyJS
using ORCA
using Rsvg
using Blink
using GeneticVariation
using ArgParse
using DelimitedFiles
using Test

export
    test_parse_main,
    format_reader,
    load_vcf,
    clean_column1!,
    load_siglist,
    #sig_list_vcf_filter,
    #genomic_range_vcf_filter,
    #load_sort_phenotype_matrix,
    #reorder_columns,
    select_columns,
    #genotype_heatmap2,
    genotype_heatmap2_new_legend,
    dp_heatmap2,
    avg_sample_dp_scatter,
    avg_variant_dp_line_chart,
    read_depth_threshhold,
    list_variant_positions_low_dp,
    list_sample_names_low_dp,
    avg_dp_variant,
    avg_dp_samples,
    jupyter_main,
    save_numerical_array,
    io_pass_filter,
    io_sig_list_vcf_filter,
    io_genomic_range_vcf_filter,
    generate_genotype_array,
    define_geno_dict,
    translate_genotype_to_num_array,
    translate_readdepth_strings_to_num_array,
    genotype_heatmap_with_groups,
    dp_heatmap2_with_groups,
    returnXY_column1!,
    pass_genomic_range_siglist_filter,
    pass_genomic_range_filter,
    pass_siglist_filter,
    genomic_range_siglist_filter,
    get_sample_names,
    sortcols_by_phenotype_matrix,
    find_group_label_indices,
    checkfor_outputdirectory,
    combined_all_genotype_array_functions,
    combined_all_read_depth_array_functions,
    combined_all_read_depth_array_functions_for_avg_dp,
    generate_chromosome_positions_for_hover_labels,
    clean_column1_chr,
    clean_column1_siglist!,
    process_plot_inputs,
    process_plot_inputs_for_grouped_data,
    returnXY_column1_siglist!,
    chromosome_label_generator,
    add_pheno_matrix_to_gt_data_for_plotting,
    add_pheno_matrix_to_dp_data_for_plotting,
    generate_hover_text_array,
    generate_hover_text_array_grouped,
    save_graphic,
    build_set_from_list,
    generate_legend_increments_ungrouped,
    generate_legend_increments_grouped


include("vcf_utils_complete.jl")
include("plot_utils.jl")
include("new_notebook_utils.jl")
include("init.jl")

end # module
