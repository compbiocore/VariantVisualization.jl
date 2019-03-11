# Plotting

##General notes: visualization options

VIVA orders all variants by chromosomal location for plotting. 

VIVA graphics are generated with PlotlyJS.jl. Graphics can be saved in *HTML*, *PDF*, *SVG*, *PNG*, and *EPS* formats. 

To create *interactive visualization* files, save VIVA's graphics in HTML format. These files are sharable and support cursor hoverlabels, zooming, panning, and PNG screen capture. Cursor hoverlabel displays genomic position, sample id, and data value for each data point in heatmap and scatter plot visualizations. We recommend saving graphics to HTML for data exploration. 

To create *publication quality, scalable graphics* for presentations and publications, we recommend saving graphics as PDF. 

## Genotype and read depth heatmaps

Plot a categorical heatmap of genotype values and a continuous value heatmap of read depth (coverage) values.

flags: `--heatmap`,`-m`

arguments: `genotype`, `read_depth`, or `genotype,read_depth`

default: `genotype,read_depth` (plots both)

```
julia VIVA -f example.vcf -m genotype
```

## Average read depth scatter plots

Generate scatter plots of average read depths across either samples or variants. Caps outlier read depth values at 100 to optimize resolution of visualization of values under 50. 


flags: `--avg_dp`

arguments: `samples`, `variants`, or `samples,variants`

default: `none`

```
julia VIVA -f example.vcf --avg_dp variants
```

## Save file format

Specify file format you wish to save all graphics as (eg. pdf, html, png). [REQUIRED]

flag: `--save_format`, `-s`  

arguments: `html`, `pdf`, `svg`, `png`, `eps`

default: `html`

```
julia VIVA -f example.vcf --avg_dp variants
``` 

## Output directory

Specify output directory for saving all graphics. If directory doesn't exist, it creates the directory within the working directory. Defaults to "output."

Select directory to save output files. If path doesn't exist, creates new directory. 

flag: `--output_directory`, `-o`

arguments: filepath

default: `output`

```
julia VIVA -f example.vcf -o my_output_directory
``` 

## Title

Specify title to display on heatmap and use as filename for saving heatmap files. Use underscores instead of spaces. Underscores will be replaced with spaces in the heatmap title.

flag: `--heatmap_title`, `-t`

arguments: title_text

default: original vcf filename

```
julia VIVA -f example.vcf -t your_heatmap_title
``` 

## Y-axis label options

Choose an option for displaying y-axis ticklabels for genomic position of variants on heatmaps and scatter plots. 
`chromosomes` separates chromosomes by adding chromosome label on the first variant of each new chromosome. 
`positions` labels every variant position (recommended only for visualizing a few variants e.g. <20)
`hoverlabels_only` no genomic position labels 

flags: `--y_axis_labels`, `-y`

arguments: `chromosomes`, `positions`, `hoverlabels_only`

default: `chromosomes`

```
julia VIVA -f example.vcf `-y` `hoverlabels_only`
```

## X-axis label options


flags: `--x_axis_labels`.`x`

arguments: if `true`, displays samples names labels on x-axis. if `false`, does not display x-axis sample labels.

default: `true`

```
julia VIVA -f example.vcf `-x`
```

## Export heatmap data as numerical array

Save input array to heatmap function with column and row labels. 
Specifically, saves numerical array of genotype or read depth values for selected variants and samples as a .csv table with genomic positions and sample names for row names and column names respectively.

flags: `--num_array`, `-n`

arguments: none, this is a positional argument.

```
julia VIVA -f example.vcf `-n` 
```
  
   
