#heatmap plots for grouped and ungrouped genotype and read depth viz

"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,filename,sample_names,gt_chromosome_labels,y_axis_label_option,x_axis_label_option,save_ext,chrom_label_info,number_rows)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option) #chr_pos_tuple_list_rev is rev because heatmap in plotly mirrors list for some reason.

            chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)
            title_no_underscores=replace(title, "_", ' ')

               trace=heatmap(
                   z = input, x = 1:size(input, 2),y = 1:size(input, 1),

                   zauto=false,zmax=3,zmin=0,

                   transpose=true,
                   colorscale =    [[0, "rgb(255,255,255)"],
                                   [0.33, "rgb(51,106,145)"],
                                   [0.66, "rgb(65,165,137)"],
                                   [1, "rgb(251,231,65)"]],

                   gridcolor = "#E2E2E2",
                   showscale = true,
                   colorbar = attr(tickvals = [0, 1, 2, 3],
                   ticktext = ["No Call (0)", "Homozygous Reference (1)", "Heterozygous Variant (2)", "Homozygous Variant (3)"])
                   );

                if y_axis_label_option == "chromosomes"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                                   ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

                  data = (trace)
                      plot(data,layout)

               elseif y_axis_label_option == "positions"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

               data = (trace)
                   plot(data,layout)

               elseif y_axis_label_option == "hover_positions"

                   layout = Layout(
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
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
   genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)

sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size = process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev)

title_no_underscores=replace(title, "_", ' ')

                  trace=heatmap(
                       z = input, x=1:size(input, 2),y=1:size(input, 1),

                       zauto=false,zmax=3,zmin=-2,

                      transpose=true,
                      colorscale = [
                                    [0, "rgb(174, 145, 255)"],
                                    [0.2, "rgb(255, 220, 145)"],
                                    [0.4, "rgb(56,25,90)"],
                                    [0.6, "rgb(51,106,145)"],
                                    [0.8, "rgb(65,165,137)"],
                                    [1, "rgb(251,231,65)"]
                                    ],

                      gridcolor = "#E2E2E2",
                      showscale = true,
                      colorbar = attr(tickvals = [-2, -1, 0, 1, 2, 3],
                      ticktext = ["Trait 2", "Trait 1", "No Call (0)", "Homozygous Reference (1)", "Heterozygous Variant (2)", "Homozygous Variant (3)"])
                      );
                      shapes = [vline(group_dividing_line)]

      if y_axis_label_option == "chromosomes"

          layout = Layout(
                          title = "$title_no_underscores",

                          xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
                        ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
    dp_heatmap2(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String}, sample_names,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of read depth data.
"""
function dp_heatmap2(input, title, chrom_label_info, sample_names,chr_pos_tuple_list_rev,y_axis_label_option)
#rename chrom labels to y_labels and y_label indices
    chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    title_no_underscores=replace(title, "_", ' ')

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

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
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                           ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
           )

          data = (trace)
              plot(data,layout)

       elseif y_axis_label_option == "positions"

           layout = Layout(
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
           )

       data = (trace)
           plot(data,layout)

       elseif y_axis_label_option == "hover_positions"

           layout = Layout(
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
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
    dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)

    sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size = process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev)

    title_no_underscores=replace(title, "_", ' ')

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        zauto=false,zmax=100,zmin=-60,

        transpose=true,

        colorscale = [
                         [0, "rgb(174, 145, 255)"],
                         [0.125, "rgb(255, 220, 145)"],
                         [0.25, "rgb(255, 255, 255)"],
                         [0.4, "rgb(153,231,255)"],
                         [0.5, "rgb(79,146,255)"],
                         [0.5625, "rgb(43,124,255)"],
                         [1, "rgb(0,64,168)"]
                     ],

                     colorbar = attr(tickvals = [-60,-40,-20,0,20,40,60,80,99],
                     title="Depth / Trait",
                     ticktext = ["Trait 2","Trait 1","No Call","0","20","40","60","80","100+"]),
                     gridcolor = "#E2E2E2",
                     showscale = true,
                     );

        shapes = [vline(group_dividing_line)]

if y_axis_label_option == "chromosomes"

    layout = Layout(
                    title = "$title_no_underscores",

                    xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                    showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
                  ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

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
    avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)
generate line chart of average read depths of each sample.
"""
function avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)

    sample_name_indices = collect(1:1:size(sample_names,2))

    trace = scatter(;x=1:size(sample_avg_list,1), y=sample_avg_list,mode="markers")
    layout = Layout(title="Average Sample Read Depth",xaxis=attr(title="Samples",ticktext=sample_names,tickvals=sample_name_indices,tickangle=45,tickfont_size=5),yaxis=attr(title="Average Read Depth"),showticklabels=false)
    plot(trace,layout)
end

"""
           avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list)
generate line chart of average read depths of each variant.
"""
function avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},chr_pos_tuple_list)

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list,1))

    trace = scatter(;x=1:size(variant_avg_list,1), y=variant_avg_list,mode="markers")
    layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions",ticktext=chr_pos_tuple_list,tickvals=chr_pos_tuple_indices,tickangle=45,tickfont_size=5,showticklabels=false),yaxis=attr(title="Average Read Depth"))
    plot(trace,layout)
end


"""
    process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)
Prepares input for heatmap plot function for both genotype and read depth plots without --group_samples flag.
"""
function process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list_rev,1))
    #chr_pos_tuple_list_rev=chr_pos_tuple_list_rev[end:-1:1,end:-1:1]

    chr_pos_tuple_list=Array{Tuple{Any,Int64}}(0)

    for i in chr_pos_tuple_list_rev
        push!(chr_pos_tuple_list, i)
    end

    sample_name_indices = collect(1:1:size(sample_names,2))
    #sample_names=sample_names[end:-1:1,end:-1:1]

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    return chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size
end

"""
    process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev)
Prepares input for heatmap plot function for both genotype and read depth plots with --group_samples flag.
"""
function process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev)

    sample_name_indices = collect(1:1:size(id_list,2))
    #id_list=id_list[end:-1:1,end:-1:1]

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    group_dividing_line=group_label_pack[3]
    group1_label=group_label_pack[4]
    group2_label=group_label_pack[5]

    chr_pos_tuple_indices = collect(1:1:size(chr_pos_tuple_list_rev,1))
    #chr_pos_tuple_list_rev=chr_pos_tuple_list_rev[end:-1:1,end:-1:1]

    chr_pos_tuple_list=Array{Tuple{Any,Int64}}(0)

    for i in chr_pos_tuple_list_rev
        push!(chr_pos_tuple_list, i)
    end

    return sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size

end
