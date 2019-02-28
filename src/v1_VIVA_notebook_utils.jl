println("Welcome to VIVA Jupyter Notebook.")
println()
println("Loading packages:")
println()

println("1/5 ArgParse")
using ArgParse

println("2/5 VIVA")
using Pkg #remove after VIVA is deployed
Pkg.activate("/Users/George/.julia/dev/VIVA") #remove after VIVA is deployed
using VIVA

println("3/5 PlotlyJS")
using PlotlyJS

println("4/5 GeneticVariation")
using GeneticVariation

println("5/5 DelimitedFiles")
using DelimitedFiles

###VCF filename
vcf_filename = "/Users/George/projects_code/ViVa.jl/test/test_files/test_4X_191.vcf"

#Variant Filters
genomic_range = "chr4:0-5000000000000" #specify genomic range of variants to select in the format chr4:0-560000000
pass_filter = "" #if "true", selects variants with Filter Status = PASS
positions_list = "" #include path and filename to list of variant positions to selects

###Sample Grouping and Selection
group_samples = "" #include path and filename to matrix to group by
select_samples = "" #include path and filename to list of sample names

###Plotting Options
heatmap = "genotype" #genotype, read_depth, or genotype,read_depth
avg_dp = "samples,variants"
y_axis_label_option = "chromosomes"#defaults to "chromosome labels". If "positions", each variant position is labeled on y axis. If "none", no y-axis ticklabels are plotted.
x_axis_label_option = "" #if true, sample ids are labelled on x axis labels

#Saving Options
heatmap_title = ""
save_format = "html"
output_directory = ""
num_array = "" #if "true", saves labeled number array of filtered heatmap input

#Passing settings to VIVA. Do not adjust.
variant_filters = genomic_range,pass_filter,positions_list
sample_selection = group_samples,select_samples
plotting_options = heatmap,avg_dp,y_axis_label_option,x_axis_label_option
saving_options = heatmap_title,save_format,output_directory,num_array


"""
jupyter_main(vcf_filename,saving_options,variant_filters,sample_selection,plotting_options)

filters, plots visualization, and saves as figure.
utilizes all global variables set in first cell of jupyter notebook
"""
function jupyter_main(vcf_filename,saving_options,variant_filters,sample_selection,plotting_options)

    genomic_range=variant_filters[1]
    pass_filter=variant_filters[2]
    positions_list=variant_filters[3]
    group_samples=sample_selection[1]
    select_samples=sample_selection[2]
    heatmap=plotting_options[1]
    avg_dp=plotting_options[2]
    y_axis_label_option=plotting_options[3]
    x_axis_labels=plotting_options[4]
    heatmap_title=saving_options[1]
    save_format=saving_options[2]
    output_directory=saving_options[3]
    num_array=saving_options[4]

    #create vcf reader object
    println("Reading $vcf_filename")
    reader = VCF.Reader(open(vcf_filename, "r"))

    #make sample names object
    sample_names = get_sample_names(reader)

    #store save format and create output directory if it doesn't exist yet
#=
    save_ext = save_format
    VIVA.checkfor_outputdirectory(output_directory)
    output_directory=output_directory
=#

    save_ext = save_format
    output_directory=output_directory

    #store plot label options
    if occursin("true",x_axis_labels)
        x_axis_label_option = true
    else
        x_axis_label_option = false
    end

    println(y_axis_label_option)

    #check for filters and apply then show stats again

    #no filter
    if !occursin(".csv",positions_list) && !occursin("chr",genomic_range) && !occursin("true",pass_filter)
        println("no filters applied. Large vcf files will take a long time to process and heatmap visualizations will not be useful at this scale.")
        println()
        println("Loading VCF file into memory for visualization")

        sub = Array{Any}(undef,0)
        for record in reader
            push!(sub,record)
        end
    end

    #all filter combos

    #pass_filter and genomic_range and list
    if occursin("true",pass_filter) && occursin("chr",genomic_range) && occursin(".csv",positions_list)
        sig_list =  load_siglist(positions_list)
        sub = VIVA.pass_chrrange_siglist_filter(vcf_filename, sig_list, genomic_range)
        number_rows = size(sub,1)
        println("Selected $number_rows variants with Filter status: PASS, that match list of chromosome positions of interest, and are within chromosome range: $genomic_range")
    end

    #pass_filter and genomic_range
    if occursin("true",pass_filter) && occursin("chr",genomic_range)
        sub = VIVA.pass_chrrange_filter(reader, genomic_range,vcf_filename)
        number_rows = size(sub,1)
        println("Selected $number_rows variants with Filter status: PASS and are within chromosome range: $genomic_range")
    end

    #pass_filter and list
    if occursin("true",pass_filter) && occursin(".csv",positions_list)
        sig_list =  load_siglist(positions_list)
        sub = VIVA.pass_siglist_filter(vcf_filename, sig_list)
        number_rows = size(sub,1)
        println("Selected $number_rows variants with Filter status: PASS and that match list of chromosome positions of interest")
    end

    #genomic_range and list
    if occursin("chr",genomic_range) && occursin(".csv",positions_list)
        sig_list =  load_siglist(positions_list)
        sub = VIVA.chrrange_siglist_filter(vcf_filename, sig_list, genomic_range)
        number_rows = size(sub,1)
        println("Selected $number_rows variants that are within chromosome range: $genomic_range and that match list of chromosome positions of interest")
    end

    #pass_filter
    if occursin("true",pass_filter) && !occursin("chr",genomic_range) && !occursin(".csv",positions_list)
        println("Only pass filter is applied. Large vcf files with many PASS variants will take a long time to process and heatmap visualizations will lose resolution at this scale unless viewed in interactive html for zooming.")
        println()
        sub = VIVA.io_pass_filter(vcf_filename)
        number_rows = size(sub,1)
        println("Selected $number_rows variants with Filter status: PASS")
        heatmap_input = "pass_filtered"
    end

    #genomic_range
    if occursin("chr",genomic_range) && !occursin("true",pass_filter) && !occursin(".csv",positions_list)
        sub = VIVA.io_chromosome_range_vcf_filter(genomic_range,vcf_filename)
        number_rows = size(sub,1)
        println("Selected $number_rows variants within chromosome range: $genomic_range")
        heatmap_input = "range_filtered"
    end

    #list
    if occursin(".csv",positions_list) && !occursin("chr",genomic_range) && !occursin("true",pass_filter)
        sig_list =  load_siglist(positions_list)
        sub = VIVA.io_sig_list_vcf_filter(sig_list,vcf_filename)
        number_rows = size(sub,1)
        println("Selected $number_rows variants that match list of chromosome positions of interest")
        heatmap_input = "positions_filtered"
    end

    println()
    println("Finished Filtering. Total time to filter:")
    println("_______________________________________________")


    if occursin("genotype",heatmap) && !occursin("read_depth",heatmap)
        gt_num_array,gt_chromosome_labels = combined_all_genotype_array_functions(sub)

        if heatmap_title != ""
            title = "Genotype_$(heatmap_title)"
        else
            bn = Base.Filesystem.basename(vcf_filename)
            title = "Genotype_$bn"
        end

        chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(gt_chromosome_labels)

        chrom_label_info = VIVA.chromosome_label_generator(gt_chromosome_labels[:,1])

        if length(split(group_samples,",")) == 2

            if select_samples != ""

                println("Make sure that select_samples list and phenotype matrix contain the same sample names.")
                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              gt_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end

            group_trait_matrix_filename=(split(group_samples,",")[1])
            trait_to_group_by = (split(group_samples,",")[2])
            println()
            println("Grouping samples by $trait_to_group_by")
            println()

            ordered_num_array,group_label_pack,pheno,id_list,trait_labels = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, gt_num_array, sample_names)

            if occursin("true",num_array)
            save_numerical_array(ordered_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            pheno_num_array,trait_label_array,chrom_label_info=add_pheno_matrix_to_gt_data_for_plotting(ordered_num_array,pheno,trait_labels,chrom_label_info,number_rows)

            graphic = VIVA.genotype_heatmap_with_groups(pheno_num_array,title,chrom_label_info,group_label_pack,id_list,chr_pos_tuple_list,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

        else

            if select_samples != ""

                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              gt_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end

            if occursin("true",num_array)
            save_numerical_array(gt_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            graphic = VIVA.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names,chr_pos_tuple_list,y_axis_label_option,x_axis_label_option)

        end

        println("Saving genotype heatmap")

        save_graphic(graphic,output_directory,save_ext,title)

    end

    if occursin("read_depth",heatmap) && !occursin("genotype",heatmap)

        dp_num_array,dp_chromosome_labels = combined_all_read_depth_array_functions(sub)

        if heatmap_title != ""
            title = "Read_Depth_$heatmap_title"
        else
            bn = Base.Filesystem.basename(vcf_filename)
            title = "Read_Depth_$bn"
        end

        chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(dp_chromosome_labels)

        chrom_label_info = VIVA.chromosome_label_generator(dp_chromosome_labels[:,1])

        if length(split(group_samples,",")) == 2

            if select_samples != ""

                println("Make sure that select_samples list and phenotype matrix contain the same sample names.")
                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              dp_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end


            group_trait_matrix_filename=(split(group_samples,",")[1])
            trait_to_group_by = (split(group_samples,",")[2])
            println()
            println("Grouping samples by $trait_to_group_by")
            println()

            ordered_dp_num_array,group_label_pack,pheno,id_list,trait_labels = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, dp_num_array, sample_names)
            dp_num_array_limited=read_depth_threshhold(ordered_dp_num_array)

             if occursin("true",num_array)
            save_numerical_array(ordered_dp_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            pheno_num_array,trait_label_array,chrom_label_info = add_pheno_matrix_to_dp_data_for_plotting(dp_num_array_limited,pheno,trait_labels,chrom_label_info,number_rows)

            graphic = VIVA.dp_heatmap2_with_groups(pheno_num_array,title,chrom_label_info,group_label_pack,id_list,chr_pos_tuple_list,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

        else

                       if select_samples != ""

                           println("Selecting samples listed in $select_samples")

                           dp_num_array,col_selectedcolumns = select_columns(select_samples,
                                                         gt_num_array,
                                                         sample_names)

                           sample_names = col_selectedcolumns
                       end

            if occursin("true",num_array)
            save_numerical_array(dp_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            dp_num_array_limited=read_depth_threshhold(dp_num_array)

            graphic = VIVA.dp_heatmap2(dp_num_array, title, chrom_label_info, sample_names,chr_pos_tuple_list,y_axis_label_option,x_axis_label_option)
        end

        println("Saving read depth heatmap")

        save_graphic(graphic,output_directory,save_ext,title)
    end

    if occursin("read_depth",heatmap) && occursin("genotype",heatmap)

        gt_num_array,gt_chromosome_labels = combined_all_genotype_array_functions(sub)

        if heatmap_title != ""
            title = "Genotype_$(heatmap_title)"
        else
            bn = Base.Filesystem.basename(vcf_filename)
            title = "Genotype_$bn"
        end

        chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(gt_chromosome_labels)

        chrom_label_info = VIVA.chromosome_label_generator(gt_chromosome_labels[:,1])

        if length(split(group_samples,",")) == 2

            if select_samples != ""

                println("Make sure that select_samples list and phenotype matrix contain the same sample names.")
                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              gt_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end

            group_trait_matrix_filename=(split(group_samples,",")[1])
            trait_to_group_by = (split(group_samples,",")[2])
            println()
            println("Grouping samples by $trait_to_group_by")
            println()

            ordered_num_array,group_label_pack,pheno,id_list,trait_labels = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, gt_num_array, sample_names)

            if occursin("true",num_array)
            save_numerical_array(ordered_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            pheno_num_array,trait_label_array,chrom_label_info=add_pheno_matrix_to_gt_data_for_plotting(ordered_num_array,pheno,trait_labels,chrom_label_info,number_rows)

            graphic = VIVA.genotype_heatmap_with_groups(pheno_num_array,title,chrom_label_info,group_label_pack,id_list,chr_pos_tuple_list,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

        else

            if select_samples != ""

                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              gt_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end

            if occursin("true",num_array)
            save_numerical_array(gt_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            graphic = VIVA.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names,chr_pos_tuple_list,y_axis_label_option,x_axis_label_option)

        end

        println("Saving genotype heatmap")

        save_graphic(graphic,output_directory,save_ext,title)

        if heatmap_title != ""
            title = "Read_Depth_$heatmap_title"
        else
            bn = Base.Filesystem.basename(vcf_filename)
            title = "Read_Depth_$bn"
        end

        dp_num_array,dp_chromosome_labels = combined_all_read_depth_array_functions(sub)

        chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(dp_chromosome_labels)

        chrom_label_info = VIVA.chromosome_label_generator(dp_chromosome_labels[:,1])

        if length(split(group_samples,",")) == 2

            if select_samples != ""

                println("Make sure that select_samples list and phenotype matrix contain the same sample names.")
                println("Selecting samples listed in $select_samples")

                gt_num_array,col_selectedcolumns = select_columns(select_samples,
                                              dp_num_array,
                                              sample_names)

                sample_names = col_selectedcolumns
            end

            group_trait_matrix_filename=(split(group_samples,",")[1])
            trait_to_group_by = (split(group_samples,",")[2])
            println()
            println("Grouping samples by $trait_to_group_by")
            println()

            ordered_dp_num_array,group_label_pack,pheno,id_list,trait_labels = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, dp_num_array, sample_names)
            dp_num_array_limited=read_depth_threshhold(ordered_dp_num_array)

            if occursin("true",num_array)
            save_numerical_array(ordered_dp_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            pheno_num_array,trait_label_array,chrom_label_info = add_pheno_matrix_to_dp_data_for_plotting(dp_num_array_limited,pheno,trait_labels,chrom_label_info,number_rows)

            graphic = VIVA.dp_heatmap2_with_groups(pheno_num_array,title,chrom_label_info,group_label_pack,id_list,chr_pos_tuple_list,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

        else


           if select_samples != ""

               println("Selecting samples listed in $select_samples")

               dp_num_array,col_selectedcolumns = select_columns(select_samples,
                                             gt_num_array,
                                             sample_names)

               sample_names = col_selectedcolumns
           end

            if occursin("true",num_array)
            save_numerical_array(dp_num_array,sample_names,chr_pos_tuple_list,title,output_directory)
            end

            dp_num_array_limited=read_depth_threshhold(dp_num_array)

            graphic = VIVA.dp_heatmap2(dp_num_array, title, chrom_label_info, sample_names,chr_pos_tuple_list,y_axis_label_option,x_axis_label_option)
        end

        println("Saving read depth heatmap")

        save_graphic(graphic,output_directory,save_ext,title)

    end

end

jupyter_main(vcf_filename,saving_options,variant_filters,sample_selection,plotting_options)
