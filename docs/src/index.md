![VIVA Logo](assets/VIVA_logo.png)

# Getting Started

# *VIVA Command Line Tool and Jupyter Notebook*

## Description

VIVA is a user-friendly command line tool built with our VariantVisualization.jl package for exploratory analysis and generation of publication quality graphics for variant analysis projects using Variant Call Format (VCF) files.

Variant selection and plotting is all executed in a single command.

We describe each of VIVA's arguments in this documentation under the Manual page.

VIVA is available as a Jupyter Notebook utility [here](https://github.com/compbiocore/VariantVisualization.jl/tree/master). Instructions for installing Jupyter, downloading VIVA Jupyter Notebook, and using the notebook are detailed in the *Jupyter Notebook* section of this documentation.

Formatting requirements for VIVA's input files are described in the Manual and clearly named examples of all user-generated input files can be found in the `/test/test_files` directory of the `VariantVisualization.jl` repository.

## General Use

To use VIVA, we recommend creating a new directory for storing your VCF file to analyze where output files will be saved. Alternatively, users may also provide paths to the VCF file and to preferred output file locations as command line arguments.

### Command Line
VIVA's general command line argument structure is as follows:

```
    julia VIVA -f file.vcf [OPTIONS]
```

From the command line or powershell, run the VIVA command line tool script which takes arguments from the command line and parses them with ArgParse.jl.

Example:

```
    julia VIVA -f example.vcf -r chr1:20000-30000000 -s pdf -m genotype,read_depth --avg_dp samples
```

To display a complete set of help instructions while using the tool, run VIVA with the help flag (`--help`, `-h`).

```
    julia VIVA -h
```

### Default Options:

By running VIVA with only a VCF filename:

```
    julia VIVA -f file.vcf
```

Default options will be used:

`--heatmap` = `genotype,read_depth`
`--save_format` = `html`
`--output_directory` = `output`
`--heatmap_title` = `vcf_filename`
`--y_axis_labels` = `chromosomes`
`--x_axis_labels` = `true`

These default settings generate a heatmap plots of genotype and read depth values of all variants for all sample ids within a VCF file.

We recommend using variant filters with most VCF files as there is too much data to plot or evaluate visually.

Specifically, we recommend visualizing fewer than 2000 variants at a time for effective visualization. However, VIVA uses memory efficient filtering and plotting and is capable of plotting >200,000 datapoints.

### Jupyter Notebook

Use the following steps to use the VIVA Jupyter Notebook utility:

1. Install Jupyter Notebook following the [platform specific instructions](https://plot.ly/python/ipython-notebook-tutorial/)

2. Download the [VIVA Jupyter Notebook](PATH) to a working directory containing your VCF file.

3. Open the Julia REPL on the command line from any directory.

4. Run `using IJulia` and then `notebook()`

5. Navigate to the directory containing the VIVA Jupyter Notebook *VIVA.ipynb* and double click to open.

6. Follow the step-by-step instructions within the notebook to generate your figures.

## Continue reading for:

* [Variant and Sample Selection](https://compbiocore.github.io/VariantVisualization.jl/stable/filtering_vcf/)

* [Plotting Options](https://compbiocore.github.io/VariantVisualization.jl/stable/plotting/)


