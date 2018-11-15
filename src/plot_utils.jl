"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option) #chr_pos_tuple_list_rev is rev because heatmap in plotly mirrors list for some reason.


            #font_size = chrom_label_info[3]
            chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

               trace=heatmap(
                   z = input, x = 1:size(input, 2),y = 1:size(input, 1),

                   zauto=false,zmax=3,zmin=0,

                   transpose=true,
                   colorscale = [[0, "rgb(255,255,255)"],
                                   [0.33, "rgb(102,212,255)"],
                                   [0.66, "rgb(61,14,105)"],
                                   [1, "rgb(236,63,69)"]],
                   gridcolor = "#E2E2E2",
                   showscale = true,
                   colorbar = attr(tickvals = [0, 1, 2, 3],
                   ticktext = ["No Call (0)", "Homozygous Reference (1)", "Heterozygous Variant (2)", "Homozygous Variant (3)"])
                   );

                if y_axis_label_option == "chromosomes"

                   layout = Layout(
                                   title = "$title",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                                   yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                                   ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

                  data = (trace)
                      plot(data,layout)

               elseif y_axis_label_option == "positions"

                   layout = Layout(
                                   title = "$title",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                                   yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
                   )

               data = (trace)
                   plot(data,layout)

               elseif y_axis_label_option == "hover_positions"


                   layout = Layout(
                                   title = "$title",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                                   yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
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

                  trace=heatmap(
                       z = input, x=1:size(input, 2),y=1:size(input, 1),

                       zauto=false,zmax=3,zmin=0,

                      transpose=true,
                      colorscale = [[0, "rgb(255,255,255)"],
                                   [0.33, "rgb(102,212,255)"],
                                   [0.66, "rgb(61,14,105)"],
                                   [1, "rgb(236,63,69)"]],
                      gridcolor = "#E2E2E2",
                      showscale = true,
                      colorbar = attr(tickvals = [0, 1, 2, 3],
                      ticktext = ["No Call (0)", "Homozygous Reference (1)", "Heterozygous Variant (2)", "Homozygous Variant (3)"])
                      );
                      shapes = [vline(group_dividing_line)]

      if y_axis_label_option == "chromosomes"

          layout = Layout(
                          title = "$title",

                          xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                          yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                          ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true),

                          shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
                        )

            data = (trace)
            plot(data,layout)

      elseif y_axis_label_option == "positions"
          layout = Layout(
                          title = "$title",

                          xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                          yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                          ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true),

                          shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
                        )

            data = (trace)
            plot(data,layout)

      elseif y_axis_label_option == "hover_positions"

        layout = Layout(
                        title = "$title",
                        xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                        ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                        yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                        ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                        shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
        )
        data = (trace)
        plot(data,layout)

      else
          println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
      end
  end

#=
                  layout = Layout(
                                  title = "$title",

                                  xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                  ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),
                                  yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                  ticktext=chr_pos_tuple_list, size=font_size,automargin=true),

                                  shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]

                                   )
                  data = (trace)
                      plot(data,layout)
              end
              =#

function pheno_matrix_heatmap(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,pheno_matrix)

      group1_index=group_label_pack[1]
      group2_index=group_label_pack[2]
      group1_dividing_line = group_label_pack[3]
      group1_label=group_label_pack[4]
      group2_label=group_label_pack[5]

      pheno_data=pheno_matrix[2:end,2:end]
      pheno_trait_labels=pheno_matrix[2:end,1]
      pheno_trait_indices=(1:size(pheno_data,1))
      id_list=id_list[end:-1:1,end:-1:1]

      trace=heatmap(
      z=pheno_data, x=1:size(pheno_data,2), y=1:size(pheno_data,1),

      zauto=false,zmax=3,zmin=0,

      transpose=true
      );
      shapes = [vline(group1_dividing_line)]

      layout = Layout(
                      xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=[group1_index, group2_index],
                      ticktext=[group1_label, group2_label]),
                      yaxis=attr(title="Phenotype", zeroline=false, tickvals=pheno_trait_indices,
                      ticktext=pheno_trait_labels),
                      shapes=shapes, yaxis_range = [1:size(pheno_data,1)], xaxis_range = [1:size(pheno_data,2)]
                       )
       data = (trace)
       plot(data,layout)
end

"""
    genotype_heatmap_with_groups_subplot(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},sample_names)
generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups_subplot(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},sample_names)

chrom_labels = chrom_label_info[1]
returnXY_column1!(chrom_labels)
chrom_label_indices = chrom_label_info[2]
font_size = chrom_label_info[3]

group1_index=group_label_pack[1]
group2_index=group_label_pack[2]
group_dividing_line=group_label_pack[3]
group1_label=group_label_pack[4]
group2_label=group_label_pack[5]
sample_names=sample_names[end:-1:1,end:-1:1]

    trace=heatmap(
         z = input, x=1:size(input,2),y=1:size(input,1),

         zauto=false,zmax=3,zmin=0,

        transpose=true,
        colorscale = [[0, "rgb(255,255,255)"],
                     [0.33, "rgb(102,212,255)"],
                     [0.66, "rgb(61,14,105)"],
                     [1, "rgb(236,63,69)"]],
        gridcolor = "#E2E2E2",
        showscale = true,
        colorbar = attr(tickvals = [0, 1, 2, 3],
        ticktext = ["No Call (0)", "Homozygous Reference (1)", "Heterozygous Variant (2)", "Homozygous Variant (3)"]),
        xaxis="x",
        yaxis="y"
        );
        shapes = [vline(group_dividing_line)]

    layout = Layout(
                    title = "$title",
                    autosize=false,width=1000, height=500,
                    xaxis=attr(showgrid=false, zeroline=false, tickvals=[group1_index, group2_index],
                    ticktext=[group1_label, group2_label]),
                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                    ticktext=chrom_labels, size=font_size),
                    shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]

                     )
    data = (trace)
        plot(data,layout)
end

"""
    dp_heatmap2(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String}, sample_names,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of read depth data.
"""
function dp_heatmap2(input, title, chrom_label_info, sample_names,chr_pos_tuple_list_rev,y_axis_label_option)

    chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        zauto=false,zmax=3,zmin=0,

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.1, "rgb(79,146,255)"],
                     [0.2, "rgb(43,124,255)"],
                     [1, "rgb(0,64,168)"]],

        colorbar = attr(tickvals = [0,20,40,60,80,99],
        title="Read Depth",
        ticktext = ["0","20","40","60","80","100+"]),

        gridcolor = "#E2E2E2",
        showscale = true,
        );

        if y_axis_label_option == "chromosomes"

           layout = Layout(
                           title = "$title",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                           yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                           ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true)
           )

          data = (trace)
              plot(data,layout)

       elseif y_axis_label_option == "positions"

           layout = Layout(
                           title = "$title",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                           yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true)
           )

       data = (trace)
           plot(data,layout)

       elseif y_axis_label_option == "hover_positions"

           layout = Layout(
                           title = "$title",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45,showticklabels=false),
                           yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false)
           )


        else
            println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
        end
    end
#=
"""
    dp_heatmap2(input, title, chrom_label_info, sample_names)
generate heatmap of read depth data.
"""
function dp_heatmap2(input, title, chrom_label_info, sample_names)

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    sample_name_indices = collect(1:1:size(sample_names,2))
    sample_names=sample_names[end:-1:1,end:-1:1]

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.1, "rgb(79,146,255)"],
                     [0.2, "rgb(43,124,255)"],
                     [1, "rgb(0,64,168)"]],
        colorbar = attr(title="Read Depth"),

        gridcolor = "#E2E2E2",
        showscale = true,
        );

    layout = Layout(
                    title = "$title",
                    xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices, ticktext=sample_names,tickfont_size=5, tickangle=45,showticklabels=false), #,showticklabels=false
                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                    ticktext=chrom_labels,size=font_size)
                    )

    data = (trace)
    plot(data,layout)
end

=#
"""
    dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev,y_axis_label_option)

    sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size = process_plot_inputs_for_grouped_data(chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},id_list,chr_pos_tuple_list_rev)

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        zauto=false,zmax=3,zmin=0,

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.1, "rgb(79,146,255)"],
                     [0.2, "rgb(43,124,255)"],
                     [1, "rgb(0,64,168)"]],
                     colorbar = attr(tickvals = [0,20,40,60,80,99],
                     title="Read Depth",
                     ticktext = ["0","20","40","60","80","100+"]),
                     gridcolor = "#E2E2E2",
                     showscale = true,
                     );
        shapes = [vline(group_dividing_line)]

if y_axis_label_option == "chromosomes"

    layout = Layout(
                    title = "$title",

                    xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                    ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true),

                    shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
                  )

      data = (trace)
      plot(data,layout)

elseif y_axis_label_option == "positions"
    layout = Layout(
                    title = "$title",

                    xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                    ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true),

                    shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
                  )

      data = (trace)
      plot(data,layout)

elseif y_axis_label_option == "hover_positions"

  layout = Layout(
                  title = "$title",
                  xaxis=attr(title="Sample ID ($(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                  ticktext=id_list, tickfont_size=5, tickangle=45,showticklabels=false),

                  yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                  ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                  shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
  )
  data = (trace)
  plot(data,layout)

else
    println("--y_axis_labels is not a valid option. Choose positions or chromosomes")
end

end

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

    trace = scatter(;x=1:size(variant_avg_list,1), y=variant_avg_list,mode="lines")
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

    chr_pos_tuple_list=Array{Tuple{Int64,Int64}}(0)

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

    chr_pos_tuple_list=Array{Tuple{Int64,Int64}}(0)

    for i in chr_pos_tuple_list_rev
        push!(chr_pos_tuple_list, i)
    end

    return sample_name_indices,id_list,chrom_labels,chrom_label_indices,font_size,group_dividing_line,group1_label,group2_label,chr_pos_tuple_indices,chr_pos_tuple_list,font_size

end
