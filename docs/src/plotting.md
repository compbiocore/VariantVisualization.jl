# Plotting

##General Notes: visualization options

Here we describe VIVA options for plotting. All plots can be generated in a single command.

VIVA orders all variants by chromosomal location for plotting.

VIVA graphics are generated with PlotlyJS.jl. Graphics can be saved in *HTML*, *PDF*, *SVG*, *PNG*, and *EPS* formats.

To create *interactive visualization* files, save VIVA's graphics in HTML format. These files are sharable and support cursor hoverlabels, zooming, panning, and PNG screen capture. Cursor hoverlabel displays genomic position, sample id, and data value for each data point in heatmap and scatter plot visualizations. We recommend saving graphics to HTML for data exploration.

To create *publication quality, scalable graphics* for presentations and publications, we recommend saving graphics as PDF.

## Genotype and read depth heatmaps

Plot a categorical heatmap of genotype values and a continuous value heatmap of read depth (coverage) values.

*flags*: `--heatmap`,`-m`

*arguments*: `genotype`, `read_depth`, or `genotype,read_depth`

default: `genotype,read_depth` (plots both)

```
julia viva -f example.vcf -m genotype
```

## Average read depth scatter plots

Generate scatter plots of average read depths across either samples or variants. Caps outlier read depth values at 100 to optimize resolution of visualization of values under 50.


*flags*: `--avg_dp`

*arguments*: `samples`, `variants`, or `samples,variants`

default: `none`

```
julia viva -f example.vcf --avg_dp variants
```

## Save file format

Specify file format you wish to save all graphics as (eg. pdf, html, png). [REQUIRED]

*flags*: `--save_format`, `-s`  

*arguments*: `html`, `pdf`, `svg`, `png`, `eps`

default: `html`

```
julia viva -f example.vcf --avg_dp variants
```

## Output directory

Specify output directory for saving all graphics. If directory doesn't exist, it creates the directory within the working directory. Defaults to "output."

Select directory to save output files. If path doesn't exist, creates new directory.

*flags*: `--output_directory`, `-o`

*arguments*: filepath

default: `output`

```
julia viva -f example.vcf -o my_output_directory
```

## Title

Specify title to display on heatmap and use as filename for saving heatmap files. Use underscores instead of spaces. Underscores will be replaced with spaces in the heatmap title.

*flags*: `--heatmap_title`, `-t`

*arguments*: title_text

default: original vcf filename

```
julia viva -f example.vcf -t your_heatmap_title
```

## Y-axis label options

Choose an option for displaying y-axis ticklabels showing the genomic position of variants on heatmaps and scatter plots.

*flags*: `--y_axis_labels`, `-y`

*arguments*: `chromosomes`, `positions`, `hoverlabels_only`

`chromosomes` separates chromosomes by adding chromosome label on the first variant of each new chromosome.
`positions` labels every variant position (recommended only for visualizing a few variants e.g. <20)
`hoverlabels_only` no genomic position labels

default: `chromosomes`

*Note*: We don't recommend using the `positions` option when visualizing samples grouped with a metadata matrix. This will show labels that are meant to be hidden that exist as an artifact of constructing the metadata trait colorbars which are sized dynamically to make up 1/20th of the plot height. If you must use the `positions` option in this scenario, we recommend editing the final plot in a program like Powerpoint to "cover up" the multitude of tick labels that will appear beside metadata trait rows. 

```
julia viva -f example.vcf `-y` `hoverlabels_only`
```

## X-axis label options

Choose an option for displaying x-axis ticklabels showing the sample id of samples included heatmaps and scatter plots.

*flags*: `--x_axis_labels`.`x`

*arguments*: if `true`, displays samples names labels on x-axis. if `false`, does not display x-axis sample labels.

default: `true`

```
julia viva -f example.vcf `-x`
```

## Export heatmap data as numerical array

Save input array to heatmap function with column and row labels.

Specifically, saves numerical array of genotype or read depth values for selected variants and samples as a .csv table with genomic positions and sample names for row names and column names respectively.

*flags*: `--num_array`, `-n`

*arguments*: none, this is a positional argument.

```
julia viva -f example.vcf `-n`
```

```
julia viva -f example.vcf `-x`
```

## Make HTML plots shareable

Make HTML plots shareable by saving HTML supporting files remotely on the internet rather than on the users local filesystem.

This allows sending HTML output files to others who don't have VIVA installed on their computers. However, in order to open HTML files saved with this option engaged, users will need to be connected to the internet.

Without calling this feature, HTML outputs can be opened without internet connection but only on the computer where the plots were generated.

*flags*: `--save_remotely`

*arguments*: none, this is a positional argument.

```
julia viva -f example.vcf `--save_remotely`
```
