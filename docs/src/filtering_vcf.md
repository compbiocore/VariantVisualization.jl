# Filtering your VCF file: Variant Record and Sample Selection

## General Notes: Extracting and Reshaping VCF Data

VIVA supports flexible filters for selecting variant records for visualization.

Additionally, the tool supports selecting and grouping samples by common traits for visualization.

Grouping samples is particularly useful for exploring phenotypic and genotypic associations, displaying differential distribution of variants between groups of samples, and identifying batch effect on coverage between groups of samples in variant analysis experiments.

## Choose a VCF file to Visualize *REQUIRED*

Specify filename of VCF file.

*flags*: `--vcf_file`, `-f`

*arguments*: Provide VCF filename (or filepath to VCF file if the file is not in the curret working directory).

*Note*: This is the *only required argument* for VIVA. If you run with none of the other options, default options will be used. These default options are described in detail below.

```
viva -f example.vcf [OPTIONS]
```

## Selecting Variant Records

VIVA offers three filters for selecting variant records to visualize from VCF files.

It is recommended to use one or a combination of these filters to reduce the number of variant records extracted from the VCF for plotting. This is recommended for reasons related to technical limitations and practical visual interpretation. The number of variant records able to be plotted is limited by both the user's available computing resources as well as the number of pixels in their display for displaying data points. While it is possible to visualize many thousands of variant records at one time with VIVA, **we recommend visualizing fewer than 2000 variants** so that all data points can be displayed that your computing resources are not overburdened. However, VIVA is capable of extracting and plotting hundreds of thousands of data points from VCF files.

### Genomic range

Select rows within a given genomic range.

*flags*: `--genomic_range`, `-r`

*arguments*: Specify genomic range within a single chromosome in format `chr4:20000000-30000000`

*Note*: To visualize genomic ranges within multiple chromosomes, you may create a batch script to run VIVA multiple times using different genomic ranges.

```
viva -f example.vcf -r chr1:20000-30000000
```

### Variant list

Select variants matching list of chromosomal positions.

*flags*: `--positions_list`, `-l`

*arguments*: Provide filename of text file formatted with two columns in .csv format as an argument. There should be a header row with "chr" and "start" in row 1 of column 1 and 2 respectively. Column 1 should contain chromosome number in the format "chr1" or "1" and should match the syntax of the VCF file (that is, if the VCF file lists chromosome numbers in the form "chrX", use "chrX" in your positions list, not "X") You can find an example of this file [here]("[here]("https://github.com/compbiocore/VariantVisualization.jl/tree/master/test/test_files/positions_list.csv")")

```
viva -f example.vcf -l "example_positions_list.txt"
```

### Pass filter

Select rows that passed filters originally set during variant calling and VCF file generation. Selects records with "PASS" in the FILTER column of the VCF file. This filter alone is often not stringent enough to reduce the number of variants for plotting and visual interpretation. For analyzing large VCF files with many "passed" filter records, use genomic range,

*flags*: `--pass_filter`, `-p`

*arguments*: This flag is a positional argument and does not take options.

```
viva -f example.vcf -p
```

## Selecting and Grouping Samples

### Group samples by sample metadata traits

Group sample columns using your sample metadata and visualize metadata attributes in a colorbar above heatmap visualizations.

*flags*: `--group_samples`, `-g`

*arguments*: This flag takes two arguments. First provide filename of sample_metadata_matrix.csv file. Second, enter the trait to group by as it appears in column 1 of the matrix.

#####  Use cases for sample grouping

VIVA supports grouping samples for visualization using any user-supplied binary metadata attributes and visualizes these in a colorbar above heatmap visualizations.

This is broadly useful for purposes such as exploring and presenting **phenotypic and genotypic associations**, **identifying batch effect on coverage**, or **visualizing differential variant incidence** between two groups of samples (such as cases and controls).

##### Input file formatting

Sample metadata matrix is a user generated input file and should be formatted in a table of sample ids and binary metadata traits (such as case,control or treatment1,treatment2 or seq_site_1,seq_site_2). An example of formatting for this table can be found [here]
("https://github.com/compbiocore/VariantVisualization.jl/tree/master/test/test_files/sample_metadata_matrix.csv").

Sample ids must match those found in the VCF file but do not need to be in the same order as they appear in the VCF header. Additionally, if the user would like to use the --select_samples option, sample ids must match the sample selection list.

Metadata traits are stored as rownames in the first column of the table and should be binary traits seperated by a comma (like "case,control"). For each metadata trait, samples should be labeled with "1" or "2" to correspond with the first and second group of the trait respectively (e.g. 1 = case and 2 = control).

This matrix should be saved as a comma delimited .csv file. Microsoft Excel is commonly used for this purpose, but sometimes creates extra delimiter characters in the output file that produce an error in VIVA. You can check to make sure the .csv file was saved properly by opening the file with a text editor such as BBEdit to inspect for and delete empty values or extra delimiter characters at the end of each row.

```
viva -f example.vcf -g sample_metadata_matrix.csv case,control
```

### Select samples to include in visualization

Select specific samples to be extracted from the VCF for visualization.

*flags*: `--select_samples`

*arguments*: Provide filename or filepath to tab delimited list of sample names to include in visualization as an argument. An example of this list can be found [here]("https://github.com/compbiocore/VariantVisualization.jl/tree/master/test/test_files/select_samples_list.txt").

*Note*: To use the sample selection feature in combination with the sample grouping feature, the sample metadata matrix must only contain the sample ids to be selected.

```
viva -f example.vcf --select_samples select_samples_list.txt
```
