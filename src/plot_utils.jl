#4) plotting functions

#a) define plotlyJS function for genotype heatmap
#call these constant variables in from module - how?
g_white = 400 #homo reference 0/0
g_red = 800 #homo variant 1/1 1/2 2/2 1/3 2/3 3/3 4/4 5/5 6/6 etc
g_pink = 600 #hetero variant 0/1 1/0 0/2 2/0 etc
g_blue = 0 #no call ./.

#red [1, "rgb(255,9,9)"]],
# yellow [1, "rgb(255,244,32)"]],

function genotype_heatmap2(x,title) #when x = array_for_plotly

    trace=heatmap(
        z = x,
        transpose=true,
        colorscale = "readdepth_colors",
        colorscale = [[0, "rgb(255,255,255)"], #choose colors and run all 6 graphics in am - replace in presentation
                     [0.5, "rgb(160,227,250)"],
                     [0.75, "rgb(52,36,242)"],
                     [1, "rgb(255,9,249)"]],
        gridcolor = "#E2E2E2",
        showscale = true
        );

    layout = Layout(
                    title = "$title",#defines title of plot
                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                    yaxis=attr(title="Chromosomal Location", zeroline=false),
                    hovermode = false
    )

    data = (trace)
        plot(data,layout) #call plot type and layout with all attributes to plot function
end

#b) define plotlyJS function for read depth heatmap

function dp_heatmap2(x, title) #when x = array_for_plotly

    #max_val=findmax(x)
    #println(max_val)

    trace=heatmap(
        z = x,

        transpose=true,
        colorscale = [[0, "rgb(153,231,255)"],
                     [0.025, "rgb(79,146,255)"],
                     [0.05, "rgb(43,124,255)"],
                     #[0.2, "rgb(0,56,147)"],
                     [1, "rgb(0,64,168)"]], #"YIGnBu","rgb(0,56,147)"

        gridcolor = "#E2E2E2",
        showscale = true
        );

    layout = Layout(
                    title = "$title",#defines title of plot
                    xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                    yaxis=attr(title="Chromosomal Location", zeroline=false,
                    hovermode = false)
    )

    data = (trace)
    plot(data,layout) #call plot type and layout with all attributes to plot function
end
