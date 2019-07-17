#heatmap plots for grouped and ungrouped genotype and read depth viz

"""
    genotype_heatmap2_new_legend(input::Array{Any,2},title::AbstractString,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option)
generate heatmap of genotype data.
"""
function genotype_heatmap2_new_legend(input,title,chrom_label_info,sample_names,chr_pos_tuple_list_rev,y_axis_label_option,x_axis_label_option) #chr_pos_tuple_list_rev is rev because heatmap in plotly mirrors list for some reason.

            chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)
            title_no_underscores=replace(title, "_"=>' ')

             hover_text_array=generate_hover_text_array(chr_pos_tuple_list,sample_names,input,"GT")

             position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,position_10,text_point_increment,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space,extra_right_margin_space_dp = generate_legend_increments_ungrouped(input)

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
                   showscale = false, #we define a custom legend using shapes and a text scatter trace
                   colorbar = attr(tickvals = [0, 1, 2, 3],
                   ticktext = ["No Call (0)", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
                   );

                   #we define a custom legend using shapes and a text scatter trace
                   trace2 = scatter(
                   ;x=[x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor+extra_right_margin_space], y=y=[position_1+text_point_increment,position_3+text_point_increment,position_5+text_point_increment,position_7+text_point_increment,position_9+text_point_increment],
                     mode="text", name="legend label",
                     textposition="middle right",
                     text=["No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant",""],
                     marker_size=300, textfont_family="Raleway, sans-serif")

                   shapes = [
                   rect(x0=legend_xaxis_anchor, y0=position_7, x1=legend_xaxis_anchor_2, y1=position_8, opacity=1.0, fillcolor="rgb(251,231,65)", line_color="black"),
                   rect(x0=legend_xaxis_anchor, y0=position_5, x1=legend_xaxis_anchor_2, y1=position_6, opacity=1.0, fillcolor="rgb(65,165,137)", line_color="black"),
                   rect(x0=legend_xaxis_anchor, y0= position_3, x1=legend_xaxis_anchor_2, y1=position_4, opacity=1.0, fillcolor="rgb(51,106,145)", line_color="black"),
                   rect(x0=legend_xaxis_anchor, y0= position_1, x1=legend_xaxis_anchor_2, y1=position_2, opacity=1.0, fillcolor="rgb(255,255,255)", line_color="black")]

                   data = [trace, trace2]

                if y_axis_label_option == "chromosomes"

                   layout = Layout(
                                   shapes=shapes,

                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickangle=45,showticklabels=x_axis_label_option),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                                   ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false)
                   )

                      plot(data,layout)

               elseif y_axis_label_option == "positions"

                   layout = Layout(shapes=shapes,
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickangle=45, showticklabels=x_axis_label_option),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false)
                   )


                   plot(data,layout)

               elseif y_axis_label_option == "hover_positions"

                   layout = Layout(shapes=shapes,
                                   title = "$title_no_underscores",
                                   xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                                   ticktext=sample_names, tickangle=45,showticklabels=x_axis_label_option),
                                   yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                                   ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false)
                   )

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

    position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,text_point_increment,trait_positions,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space = generate_legend_increments_grouped(input)

    position_1_trait_1,position_2_trait_1,position_1_trait_2,position_2_trait_2 = trait_positions

    title_no_underscores=replace(title, "_"=>' ')

    hover_text_array=generate_hover_text_array_grouped(chr_pos_tuple_list,id_list,input,"GT",number_rows)

                  trace = heatmap(
                       z = input, x=1:size(input, 2),y=1:size(input, 1),

                       hoverinfo="text",
                       text=hover_text_array,

                       zauto=false,zmax=3,zmin=-2,

                       transpose=true,

                                    colorscale = [
                                                  [0, "rgb(208, 211, 212)"], #light grey
                                                  [0.2, "rgb(151, 154, 154)"], #dark grey
                                                  [0.4, "rgb(255,255,255)"], #white
                                                  #[0.4, "rgb(56,25,90)"], #dark blue
                                                  [0.6, "rgb(51,106,145)"], #blue
                                                  [0.8, "rgb(65,165,137)"], #green
                                                  [1, "rgb(251,231,65)"] #yellow
                                                  ],

                      gridcolor = "#E2E2E2",
                      showscale = false, #we define a custom legend using shapes and a text scatter trace
                      colorbar = attr(tickvals = [-2, -1, 0, 1, 2, 3],
                      ticktext = ["Trait 2", "Trait 1", "No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
                      );

                    #we define a custom legend using shapes and a text scatter trace
                    trace2 = scatter(
                    ;x=[x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor+extra_right_margin_space,x_axis_text_anchor,x_axis_text_anchor], y=y=[position_1+text_point_increment,position_2+text_point_increment,position_3+text_point_increment,position_4+text_point_increment,position_5+text_point_increment,position_7+text_point_increment,position_8+text_point_increment],
                      mode="text", name="legend label",
                      textposition="middle right",
                      text=["No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant","", "Trait 2","Trait 1"],
                      marker_size=300, textfont_family="Raleway, sans-serif")

                    shapes = [
                    rect(x0=legend_xaxis_anchor, y0=position_4, x1=legend_xaxis_anchor_2, y1=position_5, opacity=1.0, fillcolor="rgb(251,231,65)", line_color="black"),
                    rect(x0=legend_xaxis_anchor, y0=position_3, x1=legend_xaxis_anchor_2, y1=position_4, opacity=1.0, fillcolor="rgb(65,165,137)", line_color="black"),
                    rect(x0=legend_xaxis_anchor, y0= position_2, x1=legend_xaxis_anchor_2, y1=position_3, opacity=1.0, fillcolor="rgb(51,106,145)", line_color="black"),
                    rect(x0=legend_xaxis_anchor, y0= position_1, x1=legend_xaxis_anchor_2, y1=position_2, opacity=1.0, fillcolor="rgb(255,255,255)", line_color="black"),
                    rect(x0=legend_xaxis_anchor, y0= position_7, x1=legend_xaxis_anchor_2, y1=position_8, opacity=1.0, fillcolor="rgb(151, 154, 154)", line_color="black"),
                    rect(x0=legend_xaxis_anchor, y0= position_8, x1=legend_xaxis_anchor_2, y1=position_9, opacity=1.0, fillcolor="rgb(208, 211, 212)", line_color="black")
                    ]

                     data = [trace, trace2]

      if y_axis_label_option == "chromosomes"

          layout = Layout(
                          title = "$title_no_underscores",

                          xaxis=attr(
                          title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                          showgrid=false,
                          zeroline=false,
                          tickvals=sample_name_indices,
                          ticktext=id_list,
                          tickangle=45,
                          showticklabels=x_axis_label_option),

                          yaxis=attr(
                          title="Genomic Location",
                          zeroline=false,
                          tickvals=chrom_label_indices,
                          ticktext=chrom_labels,
                          tickfont_size=font_size,
                          hovermode=true,
                          automargin=true,showgrid=false),

                          shapes=shapes,
                          yaxis_range = [1:size(input,1)],
                          xaxis_range = [1:size(input,2)]
                        )

                        plot(data,layout)

      elseif y_axis_label_option == "positions"

          layout = Layout(
                          title = "$title_no_underscores",

                          xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                          ticktext=id_list, tickangle=45,showticklabels=x_axis_label_option),

                          yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                          ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false),

                          shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
                        )

                        plot(data,layout)

      elseif y_axis_label_option == "hover_positions"

        layout = Layout(
                        title = "$title_no_underscores",
                        xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                        ticktext=id_list, tickangle=45,showticklabels=x_axis_label_option),

                        yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                        ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                        shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
        )

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

    chr_pos_tuple_indices,chr_pos_tuple_list,sample_name_indices,sample_names,chrom_labels,chrom_label_indices,font_size = process_plot_inputs(chrom_label_info,sample_names,chr_pos_tuple_list_rev)

    title_no_underscores=replace(title, "_"=>' ')

    hover_text_array=generate_hover_text_array(chr_pos_tuple_list,sample_names,input,"DP")

    position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,position_10,text_point_increment,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space,extra_right_margin_space_dp = generate_legend_increments_ungrouped(input)

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
        showscale = false, #we define a custom legend using shapes and a text scatter trace
        );

        #we define a custom legend using shapes and a text scatter trace
        trace2 = scatter(
        ;x=[x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor+extra_right_margin_space_dp],
        y=[position_1+text_point_increment,position_2+text_point_increment,position_3+text_point_increment,position_4+text_point_increment,position_5+text_point_increment,position_6+text_point_increment],
          mode="text", name="legend label",
          textposition="middle right",
          text=["0", "10", "30", "70","100+",""],
          marker_size=300, textfont_family="Raleway, sans-serif")

        shapes = [
        rect(x0=legend_xaxis_anchor, y0=position_5, x1=legend_xaxis_anchor_2, y1=position_6, opacity=1.0, fillcolor="rgb(0,64,168)", line_color="black"),
        rect(x0=legend_xaxis_anchor, y0=position_4, x1=legend_xaxis_anchor_2, y1=position_5, opacity=1.0, fillcolor="rgb(43,124,255)", line_color="black"),
        rect(x0=legend_xaxis_anchor, y0=position_3, x1=legend_xaxis_anchor_2, y1=position_4, opacity=1.0, fillcolor="rgb(79,146,255)", line_color="black"),
        rect(x0=legend_xaxis_anchor, y0= position_2, x1=legend_xaxis_anchor_2, y1=position_3, opacity=1.0, fillcolor="rgb(153,231,255)", line_color="black"),
        rect(x0=legend_xaxis_anchor, y0= position_1, x1=legend_xaxis_anchor_2, y1=position_2, opacity=1.0, fillcolor="rgb(255,255,255)", line_color="black")
        ]

        if y_axis_label_option == "chromosomes"

           layout = Layout(
                           shapes=shapes,
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chrom_label_indices,
                           ticktext=chrom_labels,tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false)
           )

          data = [trace, trace2]
              plot(data,layout)

       elseif y_axis_label_option == "positions"

           layout = Layout(
                           shapes=shapes,
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false)
           )

      data = [trace, trace2]
           plot(data,layout)

       elseif y_axis_label_option == "hover_positions"

           layout = Layout(
                           shapes=shapes,
                           title = "$title_no_underscores",
                           xaxis=attr(title="Sample ID", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickangle=45,showticklabels=x_axis_label_option),
                           yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                           ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false)
           )
          data = [trace, trace2]
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

    position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,text_point_increment,trait_positions,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space = generate_legend_increments_grouped(input)

    position_1_trait_1,position_2_trait_1,position_1_trait_2,position_2_trait_2 = trait_positions

    trace=heatmap(

        z = input, x=1:size(input, 2),y=1:size(input, 1),

        hoverinfo="text",
        text=hover_text_array,

        zauto=false,zmax=100,zmin=-60,

        transpose=true,

                     colorscale = [
                                      [0, "rgb(208, 211, 212)"], #light grey
                                      [0.125, "rgb(151, 154, 154)"], #dark grey
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
                     showscale = false, #we define a custom legend using shapes and a text scatter trace
                     );

         #we define a custom legend using shapes and a text scatter trace
         trace2 = scatter(
         ;x=[x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor,x_axis_text_anchor+extra_right_margin_space,x_axis_text_anchor,x_axis_text_anchor],
         y=[position_1+text_point_increment,position_2+text_point_increment,position_3+text_point_increment,position_4+text_point_increment,position_5+text_point_increment,position_6 + text_point_increment, position_7+text_point_increment,position_8+text_point_increment],
           mode="text", name="legend label",
           textposition="middle right",
           text=["No Call", "0", "10", "20","100+","","Trait 2","Trait 1"],
           marker_size=300, textfont_family="Raleway, sans-serif")

         data = [trace, trace2]

         shapes = [
         rect(x0=legend_xaxis_anchor, y0= position_7, x1=legend_xaxis_anchor_2, y1=position_8, opacity=1.0, fillcolor="rgb(151, 154, 154)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0= position_8, x1=legend_xaxis_anchor_2, y1=position_9, opacity=1.0, fillcolor="rgb(208, 211, 212)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0=position_5, x1=legend_xaxis_anchor_2, y1=position_6, opacity=1.0, fillcolor="rgb(37,114,242)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0=position_4, x1=legend_xaxis_anchor_2, y1=position_5, opacity=1.0, fillcolor="rgb(67,138,254)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0=position_3, x1=legend_xaxis_anchor_2, y1=position_4, opacity=1.0, fillcolor="rgb(110,182,255)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0= position_2, x1=legend_xaxis_anchor_2, y1=position_3, opacity=1.0, fillcolor="rgb(153,231,255)", line_color="black"),
         rect(x0=legend_xaxis_anchor, y0= position_1, x1=legend_xaxis_anchor_2, y1=position_2, opacity=1.0, fillcolor="rgb(255,255,255)", line_color="black")
         ]

if y_axis_label_option == "chromosomes"

    layout = Layout(
                    shapes=shapes,
                    title = "$title_no_underscores",

                    xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                    showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickangle=45,showticklabels=x_axis_label_option),

                    yaxis=attr(title="Genomic Location", zeroline=false,
                    tickvals=chrom_label_indices,ticktext=chrom_labels,
                    tickfont_size=font_size,hovermode=true,automargin=true,showgrid=false),

                    yaxis_range = [1:size(input,1)],
                    xaxis_range = [1:size(input,2)]
                  )

      plot(data,layout)

elseif y_axis_label_option == "positions"

    layout = Layout(
                    shapes=shapes,
                    title = "$title_no_underscores",

                    xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))",
                    showgrid=false, zeroline=false, tickvals=sample_name_indices,
                    ticktext=id_list, tickangle=45,showticklabels=x_axis_label_option),

                    yaxis=attr(title="Genomic Location", zeroline=false,
                    tickvals=chr_pos_tuple_indices,
                    ticktext=chr_pos_tuple_list,tickfont_size=font_size,
                    hovermode=true,automargin=true,showgrid=false),

                    yaxis_range = [1:size(input,1)],
                    xaxis_range = [1:size(input,2)]
                  )

      plot(data,layout)

elseif y_axis_label_option == "hover_positions"

  layout = Layout(
                  shapes=shapes,
                  title = "$title_no_underscores",
                  xaxis=attr(title="Sample ID (Grouped by $(group1_label) | $(group2_label))", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                  ticktext=id_list, tickangle=45,showticklabels=x_axis_label_option),

                  yaxis=attr(title="Genomic Location", zeroline=false, tickvals=chr_pos_tuple_indices,
                  ticktext=chr_pos_tuple_list,tickfont_size=font_size,hovermode=true,automargin=true,showticklabels=false),

                  yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)] #line bright white
  )
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

"""
     generate_legend_increments_ungrouped(input)
Dynamically generates positons for shapes that build categorical colorscale.
"""
function generate_legend_increments_ungrouped(input)

    number_rows = size(input,1)
    increment = number_rows/30

    number_samples = size(input,2)
    distance_from_plot_edge = number_samples/40
    width_of_color_box = number_samples/40
    legend_xaxis_anchor = number_samples + distance_from_plot_edge
    legend_xaxis_anchor_2 = legend_xaxis_anchor + width_of_color_box
    x_axis_text_anchor = legend_xaxis_anchor_2 + number_samples/40
    extra_right_margin_space=number_samples/5
    extra_right_margin_space_dp = number_samples/10

    position_1 = number_rows/2
    position_2 = position_1 + increment
    position_3 = position_2 + increment
    position_4 = position_3 + increment
    position_5 = position_4 + increment
    position_6 = position_5 + increment
    position_7 = position_6 + increment
    position_8 = position_7 + increment
    position_9 = position_8 + increment
    position_10 = position_9 + increment
    text_point_increment = increment/2
    #trait_positions = position_1_trait_1,position_2_trait_1,position_1_trait_2,position_1_trait_2

    return position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,position_10,text_point_increment,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space,extra_right_margin_space_dp
end

"""
     generate_legend_increments_grouped(input)
Dynamically generates positons for shapes that build categorical colorscale including two color boxes for traits 1 and 2
"""
function generate_legend_increments_grouped(input)

    number_rows = size(input,1)
    increment = number_rows/30

    number_samples = size(input,2)
    distance_from_plot_edge = number_samples/40
    width_of_color_box = number_samples/40
    legend_xaxis_anchor = number_samples + distance_from_plot_edge
    legend_xaxis_anchor_2 = legend_xaxis_anchor + width_of_color_box
    x_axis_text_anchor = legend_xaxis_anchor_2 + number_samples/40
    extra_right_margin_space = number_samples/5


    position_1 = number_rows/3
    position_2 = position_1 + increment
    position_3 = position_2 + increment
    position_4 = position_3 + increment
    position_5 = position_4 + increment
    position_6 = position_5 + increment
    position_7 = position_6 + increment
    position_8 = position_7 + increment
    position_9 = position_8 + increment
    #position_10 = position_9 + increment
    text_point_increment = increment/2
    position_1_trait_1 = position_9 + increment
    position_2_trait_1 = position_1_trait_1 + increment
    position_1_trait_2 = position_2_trait_1 + increment
    position_2_trait_2 = position_1_trait_2 + increment
    trait_positions = position_1_trait_1,position_2_trait_1,position_1_trait_2,position_2_trait_2

    return position_1,position_2,position_3,position_4,position_5,position_6,position_7,position_8,position_9,text_point_increment,trait_positions,legend_xaxis_anchor,legend_xaxis_anchor_2,x_axis_text_anchor,extra_right_margin_space
end
