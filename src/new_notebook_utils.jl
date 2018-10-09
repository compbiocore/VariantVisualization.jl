"""
jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,plot_types,save_format,plot_title,plot_labels)

filters, plots visualization, and saves as figure.
utilizes all global variables set in first cell of jupyter notebook
"""
function jupyter_main_new(vcf_filename::String,field_to_visualize::String,variant_filter::Array{String,1},sample_filter::Array{String,1},plot_types::Array{String,1},save_format::String,plot_title::String,output_directory::String)

println("loading packages...")

println("finished loading packages")
#show stats
#=
number_records = nrecords(vcf_filename)
number_samples = nsamples(vcf_filename)

println("_______________________________________________")
println()
println("Summary Statistics of $(parsed_args["vcf_file"])")
println()
println("number of records: $number_records")
println("number of samples: $number_samples")
println("_______________________________________________")
println()
=#
#create vcf reader object
println("Reading $vcf_filename")
reader = VCF.Reader(open(vcf_filename, "r")) #reader = VCF.Reader(open("test_4X_191.vcf", "r"))

#make sample names object
sample_names = get_sample_names(reader)

#check for filters and apply then show stats again

if typeof(variant_filter) == nothing
    println("no filters applied. Large vcf files will take a long time to process and heatmap visualizations will not be useful at this scale.")
    println()
    println("Loading VCF file into memory for visualization")

    sub = Array{Any}(0)
    for record in reader
        push!(sub,record)
    end
end

if typeof(variant_filter) != String

    for i = 1:size(variant_filter,1)

        if variant_filter[i] == "pass_only"

            sub = ViVa.io_pass_filter(reader)
            number_rows = size(sub,1)
            println("selected $number_rows variants that with Filter status: PASS")
            heatmap_input = "pass_filtered"

        elseif variant_filter[i] == "range"

             chr_range = variant_filter[i+1]
             println(chr_range)
             sub = ViVa.io_chromosome_range_vcf_filter(chr_range,reader)
             number_rows = size(sub,1)
             println("selected $number_rows variants within chromosome range: $chr_range)")
             heatmap_input = "range_filtered"

        elseif variant_filter[i] == "list"
            siglist_file = variant_filter[i+1]
            sub = ViVa.io_sig_list_vcf_filter(siglist_file,reader)
            number_rows = size(sub,1)
            println("selected $number_rows variants that match list of chromosome positions of interest")
            heatmap_input = "positions_filtered"

        end
    end

else

    if variant_filter == "pass_only"
        sub = ViVa.io_pass_filter(reader)
        number_rows = size(sub,1)
        println("selected $number_rows variants that with Filter status: PASS")
        heatmap_input = "pass_filtered"

    elseif variant_filter == "range"
        chr_range = variant_filter[i+1]
        sub = ViVa.io_chromosome_range_vcf_filter(chr_range,reader)
        number_rows = size(sub,1)
        println("selected $number_rows variants within chromosome range: $chr_range)")
        heatmap_input = "range_filtered"

    elseif variant_filter == "list"
        siglist_file = variant_filter[i+1]
        sub = ViVa.io_sig_list_vcf_filter(siglist_file,reader)
        number_rows = size(sub,1)
        println("selected $number_rows variants that match list of chromosome positions of interest")
        heatmap_input = "positions_filtered"

    end
end


#check for samples filters and apply

if typeof(sample_filter) != String

    for i = 1:size(sample_filter,1)
        if sample_filter[i] == "reorder_columns"

            list =sample_filter[i+1]
            key = sample_filter[i+2]
            println("sorting columns by $key in $list")

        elseif sample_filter[i] == "select_columns"
            id_list = sample_filter[i+1]
            println("selecting columns to match $id_list")

        end
    end

else
    if sample_filter == "reorder_columns"
        list =sample_filter[i+1]
        key = sample_filter[i+2]
        println("sorting columns by $key in $list")

    elseif sample_filter == "select_columns"
        id_list = sample_filter[i+1]
        println("selecting columns to match $id_list")
    end
end

#check plot types and generate, save, and display image

#plots_to_generate = split(plot_types,",")
#plots_to_generate = split(plots_to_generate1,",")
#(plots_to_generate)

println(sample_filter)

if field_to_visualize == "genotype"

    #genotype_array = generate_genotype_array(sub,"GT")

              #clean_column1!(genotype_array)
              #genotype_array=ViVa.sort_genotype_array(genotype_array)

              #geno_dict = define_geno_dict()
              #gt_num_array,gt_chromosome_labels = translate_genotype_to_num_array(genotype_array, geno_dict)
              gt_num_array,gt_chromosome_labels = combined_all_genotype_array_functions(sub)

              if plot_title != nothing
                  title = plot_title
              else
                  title = "Genotype_$vcf_filename)"
              end

              chrom_label_info = ViVa.chromosome_label_generator(gt_chromosome_labels[:,1])

              #select_columns and reorder_columns
              #for i = 1:size(sample_filter,1)
              if in("reorder_columns",sample_filter)
                  list =sample_filter[2]#list =sample_filter[i+1]
                  key = sample_filter[3]#key = sample_filter[i+2]
                  println(key)
                  println("sorting columns by $key in $list")
                  gt_num_array,group_label_pack = sortcols_by_phenotype_matrix(list, key, gt_num_array, sample_names)
                  graphic = ViVa.genotype_heatmap_with_groups(gt_num_array,title,chrom_label_info,group_label_pack,sample_names)

              elseif in("select_columns",sample_filter)
                  id_list = sample_filter[2] #id_list = sample_filter[i+1]
                  println("selecting columns to match $id_list")
                  gt_num_array = select_columns(id_list, gt_num_array, sample_names)
                  graphic = ViVa.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names)

              else
                  graphic = ViVa.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names)

              end
          #end

              #graphic = ViVa.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names)
              ViVa.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names)
              PlotlyJS.savefig(graphic, joinpath("$output_directory","$title.$save_format"))
              return graphic

elseif field_to_visualize == "read_depth"

    println("plotting read depth")

    read_depth_array = ViVa.generate_genotype_array(sub,"DP")

      clean_column1!(read_depth_array)
      read_depth_array=ViVa.sort_genotype_array(read_depth_array)

      dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array(read_depth_array)

  if in("heatmap",plot_types)

      if plot_title != nothing
          title = plot_title
      else
          title = "Read_depth_$vcf_filename"
      end

      chrom_label_info = ViVa.chromosome_label_generator(dp_chromosome_labels[:,1])

      if in("reorder_columns",sample_filter)
          list =sample_filter[2]#list =sample_filter[i+1]
          key = sample_filter[3]#key = sample_filter[i+2]
          println(key)
          println("sorting columns by $key in $list")
          dp_num_array,group_label_pack = sortcols_by_phenotype_matrix(list, key, dp_num_array, sample_names)
          dp_num_array_limited=read_depth_threshhold(dp_num_array) #add option for this, default is true, option to turn off
          dp_graphic = ViVa.dp_heatmap2_with_groups(dp_num_array_limited,title,chrom_label_info,group_label_pack)

      elseif in("select_columns",sample_filter)
          id_list = sample_filter[2] #id_list = sample_filter[i+1]
          println("selecting columns to match $id_list")
          dp_num_array = select_columns(id_list, dp_num_array, sample_names)
          dp_num_array_limited=read_depth_threshhold(dp_num_array) #add option for this, default is true, option to turn off
          dp_graphic = ViVa.dp_heatmap2_with_groups(dp_num_array_limited,title,chrom_label_info,group_label_pack)

      else
          dp_num_array_limited=read_depth_threshhold(dp_num_array) #add option for this, default is true, option to turn off
          dp_graphic = ViVa.dp_heatmap2(dp_num_array_limited,title,chrom_label_info,group_label_pack)

      end

      #dp_num_array_limited=read_depth_threshhold(dp_num_array) #add option for this, default is true, option to turn off
      #graphic = ViVa.dp_heatmap2(dp_num_array_limited,title,chrom_label_info)
      ViVa.dp_heatmap2(dp_num_array_limited,title,chrom_label_info)
      PlotlyJS.savefig(dp_graphic, joinpath("$output_directory","$title.$save_format"))
      return dp_graphic

  end

  if in("variant_line_chart",plot_types)

      avg_list = ViVa.avg_dp_variant(dp_num_array)
      list = list_variant_positions_low_dp(avg_list, dp_chromosome_labels)
      writedlm(joinpath("$output_directory","Variant_positions_with_low_dp.txt"),list,",")

      #println("The following samples have read depth of less than 15: $list")
      graphic = avg_variant_dp_line_chart(avg_list)
      avg_variant_dp_line_chart(avg_list)
      PlotlyJS.savefig(graphic, joinpath("$output_directory","Average_Variant_Read_Depth.$save_format"))
      return graphic


  end

  if in("sample_line_chart",plot_types)

       avg_list = ViVa.avg_dp_samples(dp_num_array)
       list = ViVa.list_sample_names_low_dp(avg_list, sample_names)
       writedlm(joinpath("$output_directory","Samples_with_low_dp.txt"),list, ",")
       #println("The following samples have read depth of under 15: $list")
       avg_sample_dp_line_chart(avg_list)
       graphic = avg_sample_dp_line_chart(avg_list)
       PlotlyJS.savefig(graphic, joinpath("$output_directory" ,"Average Sample Read Depth.$save_format"))
       return graphic

  end

end

#return graphic_geno_heatmap,graphic_dp_heatmap,graphic_vlc,graphic_slc

close(reader)


end #end of jupyter_main()
