#heatmap plots for grouped and ungrouped genotype and read depth viz

"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,filename,sample_names,gt_chromosome_labels,y_axis_label_option,x_axis_label_option)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option) #chr_pos_tuple_list_rev is rev because heatmap in plotly mirrors list for some reason.

            chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)
            title_no_underscores=replace(title, "_"=>' ')

             hover_text_array=generate_hover_text_array(chr_pos_tuple_list,sample_names,input,"GT")

               trace=heatmap(
                   z = input, x = 1:size(input, 2),y = 1:size(input, 1),

                   hoverinfo="text",
                   text=hover_text_array,

                   zauto=false,zmax=3,zmin=0,

                   transpose=true,
                   colorscale =    [[0, "rgb(255,255,255)"],
                                   [0.33, "rgb(51,106,145)"],
                                   [0.66, "rgb(65,165,137)"],
                                   [1, "rgb(251,231,65)"]],

                   gridcolor = "#E2E2E2",
                   showscale = true,
                   colorbar = attr(tickvals = [0, 1, 2, 3],
                   ticktext = ["No Call (0)", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
                   );

                if y_axis_label_option == "chromosomes"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=true),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                                   ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

                  data = (trace)
                      plot(data,layout)

               elseif y_axis_label_option == "positions"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45, showticklabels=x_axis_label_option),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

               data = (trace)
                   plot(data,layout)

               elseif y_axis_label_option == "hover_positions"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false)
                   )

               data = (trace)
                   plot(data,layout)

                else
                    println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
                end
end

"""
   genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)
generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size = process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,trait_label_array)

title_no_underscores=replace(title, "_"=>' ')

hover_text_array=generate_hover_text_array_grouped(chr_pos_tuple_list,id_list,input,"GT",number_rows)

h_line_index_list = [number_rows + 0.5]
#=
function generate_hline_indices(number_rows,input)


    h_line_index_list = [number_rows + 0.5]

    traits=unique(trait_label_array)
    pheno_row_multiplyer=0.05*(number_rows)
    pheno_row_multiplyer=(pheno_row_multiplyer)
    pheno_matrix_size=size(input,1)-number_rows
    number_trait_rows=round(pheno_row_multiplyer/length(traits))+1
    println(number_trait_rows)

    #pheno_row_multiplyer=(0.05*(number_rows))
    #println(pheno_row_multiplyer)

    #pheno_row_multiplyer=(pheno_row_multiplyer/(size(input,1)-number_rows))
    #println(pheno_row_multiplyer)

    #number_trait_rows=((size(input,1)-number_rows)/pheno_row_multiplyer)
    #println(number_trait_rows)

    #number_trait_rows=round(number_trait_rows/(length(traits)))-1
    #println(number_trait_rows)



    for n=1:length(traits)
    println(number_trait_rows)
        trait_line_index = h_line_index_list[n]+number_trait_rows
         push!(h_line_index_list,trait_line_index)
    end

    println(h_line_index_list)

    #h_line_index_list = [1484,1570]

return h_line_index_list
end

h_line_index_list = generate_hline_indices(number_rows,input)
=#

                  trace=heatmap(
                       z = input, x=1:size(input, 2),y=1:size(input, 1),

                       hoverinfo="text",
                       text=hover_text_array,

                       zauto=false,zmax=3,zmin=-2,

                      transpose=true,

                      #=
                      colorscale = [
                      [0, "rgb(210, 213, 219)"], #light grey
                      [0.2, "rgb(123, 125, 130)"], #dark grey
                                    [0.4, "rgb(255,255,255)"], #white
                                    #[0.4, "rgb(56,25,90)"], #dark blue
                                    [0.6, "rgb(51,106,145)"], #blue
                                    [0.8, "rgb(65,165,137)"], #green
                                    [1, "rgb(251,231,65)"] #yellow
                                    ],
=#

                                    colorscale = [
                                                  [0, "rgb(174, 145, 255)"], #purple
                                                  [0.2, "rgb(255, 220, 145)"], #gold
                                                  [0.4, "rgb(255,255,255)"], #white
                                                  #[0.4, "rgb(56,25,90)"], #dark blue
                                                  [0.6, "rgb(51,106,145)"], #blue
                                                  [0.8, "rgb(65,165,137)"], #green
                                                  [1, "rgb(251,231,65)"] #yellow
                                                  ],

                      gridcolor = "#E2E2E2",
                      showscale = true,
                      colorbar = attr(tickvals = [-2, -1, 0, 1, 2, 3],
                      ticktext = ["Trait 2", "Trait 1", "No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
                      );

                      shapes_definition = [vline(group_dividing_line)]

                     # h_line_index_list = [h_line_index, h_line_index+50]

                      for i in h_line_index_list

                          horizontal = hline(i,line=attr(color="white"))
                          shapes_definition=push!(shapes_definition,horizontal)

                      end

                      shapes=shapes_definition


                #    shapes = [vline(group_dividing_line), hline(h_line_index,line=attr(color="white"))]
                #    println(shapes)

      if y_axis_label_option == "chromosomes"

          layout = Layout(
                          title = "$title_no_underscores",

                          xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                          yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                          ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true),

                          shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
                        )

            data = (trace)
            plot(data,layout)

      elseif y_axis_label_option == "positions"
          layout = Layout(
                          title = "$title_no_underscores",

                          xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                          yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                          ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true),

                          shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
                        )

            data = (trace)
            plot(data,layout)

      elseif y_axis_label_option == "hover_positions"

        layout = Layout(
                        title = "$title_no_underscores",
                        xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                        ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                        yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                        ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                        shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
        )
        data = (trace)
        plot(data,layout)

      else
          println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
      end
  end

"""
    dp_heatmap2(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String}, sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option)
generate heatmap of read depth data.
"""
function dp_heatmap2(input, title, chrom_label_info, sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option)
#rename chrom labels to y_labels and y_label indices
    chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    title_no_underscores=replace(title, "_"=>' ')

     hover_text_array=generate_hover_text_array(chr_pos_tuple_list,sample_names,input,"DP")

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

       hoverinfo="text",
       text=hover_text_array,

        zauto=false,zmax=100,zmin=-20,

        transpose=true,
        colorscale = [
                     [0, "rgb(255,255,255)"],
                     [0.2, "rgb(153,231,255)"],
                     [0.3, "rgb(79,146,255)"],
                     [0.41, "rgb(43,124,255)"],
                     [1, "rgb(0,64,168)"]
                     ],

        colorbar = attr(tickvals = [-20,0,20,40,60,80,99],
        title="Read Depth",
        ticktext = ["No Call","0","20","40","60","80","100+"]),

        gridcolor = "#E2E2E2",
        showscale = true,
        );

        if y_axis_label_option == "chromosomes"

           layout = Layout(
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                           ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
           )

          data = (trace)
              plot(data,layout)

       elseif y_axis_label_option == "positions"

           layout = Layout(
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
           )

       data = (trace)
           plot(data,layout)

       elseif y_axis_label_option == "hover_positions"

           layout = Layout(
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false)
           )
           data = (trace)
               plot(data,layout)

        else
            println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
        end
end

"""
    dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option,trait_label_array,x_axis_label_option,number_rows)

    sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size = process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,trait_label_array)

    title_no_underscores=replace(title, "_"=>' ')

     hover_text_array=generate_hover_text_array_grouped(chr_pos_tuple_list,id_list,input,"DP",number_rows)

    h_line_index = (number_rows + 0.5)

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        hoverinfo="text",
        text=hover_text_array,

        zauto=false,zmax=100,zmin=-60,

        transpose=true,

#=

        colorscale = [
        [0, "rgb(210, 213, 219)"], #light grey
        [1.25, "rgb(123, 125, 130)"], #dark grey
                         [0.25, "rgb(255, 255, 255)"],
                         [0.4, "rgb(153,231,255)"],
                         [0.5, "rgb(79,146,255)"],
                         [0.5625, "rgb(43,124,255)"],
                         [1, "rgb(0,64,168)"]
                     ],

                =#
                     colorscale = [
                                      [0, "rgb(174, 145, 255)"], #purple
                                      [0.125, "rgb(255, 220, 145)"], #gold
                                      [0.25, "rgb(255, 255, 255)"],
                                      [0.4, "rgb(153,231,255)"],
                                      [0.5, "rgb(79,146,255)"],
                                      [0.5625, "rgb(43,124,255)"],
                                      [1, "rgb(0,64,168)"]
                                  ],


                     colorbar = attr(tickvals = [-60, -40, -20, 0, 20, 40, 60, 80, 99],
                     title="Depth / Trait",
                     ticktext = ["Trait 2","Trait 1","No Call","0","20","40","60","80","100+"]),
                     gridcolor = "#E2E2E2",
                     showscale = true,
                     );

                     shapes = [vline(group_dividing_line), hline(h_line_index,line=attr(color="white"))]

if y_axis_label_option == "chromosomes"

    layout = Layout(
                    title = "$title_no_underscores",

                    xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                    showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                    yaxis=attr(title="Genomic Location", zeroline=false,
                    tickvals=chrom_label_indices,ticktext=chrom_labels,
                    tickfont_size=font_size,hovermode=true,automargin=true),

                    shapes=shapes, yaxis_range = [1:size(input,1)],
                    xaxis_range = [1:size(input,2)]
                  )

      data = (trace)
      plot(data,layout)

elseif y_axis_label_option == "positions"
    layout = Layout(
                    title = "$title_no_underscores",

                    xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                    showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                    yaxis=attr(title="Genomic Location", zeroline=false,
                    tickvals=chr_pos_tuple_indices,
                    ticktext=chr_pos_tuple_list,tickfont_size=font_size,
                    hovermode=true,automargin=true),

                    shapes=shapes, yaxis_range = [1:size(input,1)],
                    xaxis_range = [1:size(input,2)] #line bright white
                  )

      data = (trace)
      plot(data,layout)

elseif y_axis_label_option == "hover_positions"

  layout = Layout(
                  title = "$title_no_underscores",
                  xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                  ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=x_axis_label_option),

                  yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                  ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                  shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
  )
  data = (trace)
  plot(data,layout)

else
    println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
end

end

#Scatter plots for average read depth viz
"""
    avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names,x_axis_label_option)
generate line chart of average read depths of each sample.
"""
function avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names,x_axis_label_option)

    sample_name_indices = collect(1:1:size(sample_names,2))

    trace = scatter(;x=1:size(sample_avg_list,1), y=sample_avg_list,mode="markers")


    layout = Layout(title="Average Sample Read Depth",
    xaxis=attr(title="Samples",ticktext=sample_names,tickvals=sample_name_indices,
    tickangle=45,tickfont_size=5,showticklabels=x_axis_label_option),
    yaxis=attr(title="Average Read Depth"),showticklabels=false)

    plot(trace,layout)
end

"""
    avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list,y_axis_label_option,chrom_label_info)
generate line chart of average read depths of each variant.
"""
function avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list,y_axis_label_option,chrom_label_info)

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list,1))
    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    trace = scatter(;x=1:size(variant_avg_list,1), y=variant_avg_list,mode="markers")

    if y_axis_label_option == "positions"
        layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions",
        ticktext=chr_pos_tuple_list,
        tickvals=chr_pos_tuple_indices,
        tickangle=45,showticklabels=true),yaxis=attr(title="Average Read Depth"))
    elseif y_axis_label_option == "chromosomes"
        layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions",
        tickvals=chrom_label_indices,
        ticktext=chrom_labels,
        showticklabels=true),yaxis=attr(title="Average Read Depth"))
    elseif y_axis_label_option == "hover_positions"
        layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions",
        ticktext=chr_pos_tuple_list,
        tickvals=chr_pos_tuple_indices,
        tickangle=45,showticklabels=false,hovermode=true),yaxis=attr(title="Average Read Depth"))
    else
        println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
    end

    plot(trace,layout)

end

"""
    process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)
Prepares input for heatmap plot function for both genotype and read depth plots without --group_samples flag.
"""
function process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list_rev,1))

    #chr_pos_tuple_list=Array{Tuple{Any,Int64}}(undef,0) #use this if tuple to string doesnt work #version
    chr_pos_tuple_list=Array{String}(undef,0)

    for i in chr_pos_tuple_list_rev
        push!(chr_pos_tuple_list, i)
    end

    sample_name_indices = collect(1:1:size(sample_names,2))

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    return chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size
end

"""
    process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,trait_label_array)
Prepares input for heatmap plot function for both genotype and read depth plots with --group_samples flag.
"""
function process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,trait_label_array)

    sample_name_indices = collect(1:1:size(id_list,2))

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    group_dividing_line=group_label_pack[3]
    group1_label=group_label_pack[4]
    group2_label=group_label_pack[5]

    chr_pos_tuple_list=Array{String}(undef,0)

    for i in chr_pos_tuple_list_rev
        push!(chr_pos_tuple_list, i)
    end

    y_labels=vcat(chr_pos_tuple_list,trait_label_array)

    y_labels_indices = collect(1:1:size(y_labels,1))

    return sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,y_labels_indices,y_labels,font_size

end
