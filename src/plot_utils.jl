#=**modify colors in arguments***
sortby A then by B
fix dp parsing with threshhold
open notebook() from CLI
phenomatrix sort by name case:control
=#

#a) define plotlyJS function for genotype heatmap
#call these constant variables in from module - how?

g_white = 400 #homo reference 0/0
g_red = 800 #homo variant 1/1 1/2 2/2 1/3 2/3 3/3 4/4 5/5 6/6 etc
g_pink = 600 #hetero variant 0/1 1/0 0/2 2/0 etc
g_blue = 0 #no call ./.

#=

tick font -
tickfont=dict(
        family='Old Standard TT, serif',
        size=8,
        color='black'

        =#

"""
    genotype_heatmap2(input::Array{Any,2},title::AbstractString,chrom_label_info,sample_names)
generate heatmap of genotype data.
"""
function genotype_heatmap2(input,title,chrom_label_info,sample_names)

            chrom_labels = chrom_label_info[1]
            returnXY_column1!(chrom_labels)
            chrom_label_indices = chrom_label_info[2]
            font_size = chrom_label_info[3]
            sample_name_indices = collect(1:1:size(sample_names,2))

           trace=heatmap(
               z = input, x = 1:size(input, 2),y = 1:size(input, 1),
              # text = sample_names,
               #hoverinfo = "text",
               transpose=true,
               #colorscale = "readdepth_colors",
               colorscale = [[0, "rgb(255,255,255)"], #choose colors and run all 6 graphics in am - replace in presentation
                            [0.5, "rgb(102,212,255)"],  #"rgb(160,227,250)" original blue
                            [0.75, "rgb(61,14,105)"],
                            [1, "rgb(236,63,69)"]],   #[1, "rgb(255,9,249)"]] original pink
               gridcolor = "#E2E2E2",
               showscale = true,
               colorbar = attr(tickvals = [0, 400, 600, 800],
               ticktext = ["No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
               );

               #trace2=scatter(shapes = line(x0=100 ,y0=0, x1=100, y1=size(input,1))# color="rgb(55, 128, 191)", width=3))

           layout = Layout(
                           title = "$title",#defines title of plot
                           xaxis=attr(title="Sample Number", showgrid=false, zeroline=false, tickvals=sample_name_indices,
                           ticktext=sample_names, tickfont_size=5, tickangle=45),
                           yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                           ticktext=chrom_labels,tickfont_size=font_size,hovermode=true)

           )
           data = (trace)
               plot(data,layout) #call plot type and layout with all attributes to plot function
       end

"""
           genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},sample_names)
       generate heatmap of genotype data.
"""
function genotype_heatmap_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1},sample_names)

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    group1_index=group_label_pack[1]
    group2_index=group_label_pack[2]
    group_dividing_line=group_label_pack[3]
    group1_label=group_label_pack[4]
    group2_label=group_label_pack[5]

                  trace=heatmap(
                       z = input, x=1:size(input, 2),y=1:size(input, 1),

                      transpose=true,
                      #colorscale = "readdepth_colors",
                      colorscale = [[0, "rgb(255,255,255)"], #choose colors and run all 6 graphics in am - replace in presentation
                                   [0.5, "rgb(160,227,250)"],
                                   [0.75, "rgb(52,36,242)"],
                                   [1, "rgb(236,63,69)"]],
                      gridcolor = "#E2E2E2",
                      showscale = true,
                      colorbar = attr(tickvals = [0, 400, 600, 800],
                      ticktext = ["No Call", "Homozygous Reference", "Heterozygous Variant", "Homozygous Variant"])
                      );
                      shapes = [vline(group_dividing_line)]

                  layout = Layout(
                                  title = "$title",#defines title of plot

                                  xaxis=attr(title="Sample Number", showgrid=false, zeroline=false, tickvals=[group1_index, group2_index],
                                  ticktext=[group1_label, group2_label]),

                                  yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                                  ticktext=chrom_labels, size=font_size),
                                  shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]

                                   )
                  data = (trace)
                      plot(data,layout) #call plot type and layout with all attributes to plot function
              end


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

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.1, "rgb(79,146,255)"],
                     [0.2, "rgb(43,124,255)"],#0.05
                     #[0.2, "rgb(0,56,147)"],
                     [1, "rgb(0,64,168)"]], #"YIGnBu","rgbrgb(0,56,147)"
        colorbar = attr(title="Read Depth"),

        gridcolor = "#E2E2E2",
        showscale = true,
        );

    layout = Layout(
                    title = "$title",#defines title of plot
                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false, tickvals=sample_name_indices, ticktext=sample_names,tickfont_size=5, tickangle=45), #,showticklabels=false
                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                    ticktext=chrom_labels,size=font_size)
                    )

    data = (trace)
    plot(data,layout) #call plot type and layout with all attributes to plot function
end

"""
           dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1})
generate heatmap of read depth data with grouped samples.
"""
function dp_heatmap2_with_groups(input::Array{Int64,2},title::String,chrom_label_info::Tuple{Array{String,1},Array{Int64,1},String},group_label_pack::Array{Any,1})

    chrom_labels = chrom_label_info[1]
    returnXY_column1!(chrom_labels)
    chrom_label_indices = chrom_label_info[2]
    font_size = chrom_label_info[3]

    group1_index=group_label_pack[1]
    group2_index=group_label_pack[2]
    group_dividing_line=group_label_pack[3]
    group1_label=group_label_pack[4]
    group2_label=group_label_pack[5]

    trace=heatmap(
        z = input, x=1:size(input, 2),y=1:size(input, 1),

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.1, "rgb(79,146,255)"],
                     [0.2, "rgb(43,124,255)"],
                     #[0.2, "rgb(0,56,147)"],
                     [1, "rgb(0,64,168)"]], #"YIGnBu","rgbrgb(0,56,147)"
        colorbar = attr(title="Read Depth"),

        gridcolor = "#E2E2E2",
        showscale = true,
        );
    shapes = [vline(group_dividing_line)]

    layout = Layout(
                    title = "$title",#defines title of plot

                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false, tickvals=[group1_index, group2_index],
                    ticktext=[group1_label, group2_label]),

                    yaxis=attr(title="Chromosomal Location", zeroline=false, tickvals=chrom_label_indices,
                    ticktext=chrom_labels, size=font_size),
                    shapes=shapes, yaxis_range = [1:size(input,1)], xaxis_range = [1:size(input,2)]
                    )

    data = (trace)
    plot(data,layout)
end

"""
           avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)
generate line chart of average read depths of each sample.
"""
function avg_sample_dp_scatter(sample_avg_list::Array{Float64,1},sample_names)

    sample_name_indices = collect(1:1:size(sample_names,2))

    trace = scatter(;x=1:size(sample_avg_list,1), y=sample_avg_list,mode="markers") #, mode="lines"
    layout = Layout(title="Average Sample Read Depth",xaxis=attr(title="Samples",ticktext=sample_names,tickvals=sample_name_indices,tickangle=45,tickfont_size=5),yaxis=attr(title="Average Read Depth"),showticklabels=false)
    plot(trace,layout)
end

"""
           avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},sample_names)
generate line chart of average read depths of each variant.
"""
function avg_variant_dp_line_chart(variant_avg_list::Array{Float64,1},sample_names)

    sample_name_indices = collect(1:1:size(sample_names,2))

    trace = scatter(;x=1:size(variant_avg_list,1), y=variant_avg_list,mode="lines") #,text="test_text" , mode="lines+text"
    layout = Layout(title="Average Variant Read Depth",xaxis=attr(title="Variant Positions"),yaxis=attr(title="Average Read Depth"),tickangle=45,tickfont_size=5)
    plot(trace,layout)
end
