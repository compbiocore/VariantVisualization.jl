#4) plotting functions

#a) define plotlyJS function for genotype heatmap

function genotype_heatmap2(x,title) #when x = array_for_plotly

    trace=heatmap(
        z = x,
        transpose=true,
        #colorscale = "Picnic",
        colorscale = "Rainbow",
        gridcolor = "#E2E2E2",
        showscale = true
        );

    layout = Layout(
                    title = "$title",#defines title of plot
                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                    yaxis=attr(title="Chromosomal Location", zeroline=false)
    )

    data = (trace)
    plot(data,layout) #call plot type and layout with all attributes to plot function
end

#b) define plotlyJS function for read depth heatmap

function dp_heatmap2(x,title) #when x = array_for_plotly

    trace=heatmap(
        z = x,
        transpose=true,
        colorscale = "YIGnBu",
        gridcolor = "#E2E2E2",
        showscale = true
        );

    layout = Layout(
                    title = "$title",#defines title of plot
                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                    yaxis=attr(title="Chromosomal Location", zeroline=false)
    )

    data = (trace)
    plot(data,layout) #call plot type and layout with all attributes to plot function
end
