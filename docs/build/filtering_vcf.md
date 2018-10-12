
<a id='VCF-Record-and-Sample-Manipulation-1'></a>

# VCF Record and Sample Manipulation


```
ViVa.jl provides tools for filtering records as well as selecting and grouping samples by common traits for visualization. Users can access these tools by using ViVa.jl on the command line using the following arguments or single character flags when available for convenience.
```


<a id='Reading-a-VCF-1'></a>

# Reading a VCF


```
Specify filename of VCF file. [REQUIRED]
"--vcf_file", "-f"
... -f example.vcf ..

Print number of records and samples in VCF file.
"--show_stats"
... --show_stats ...
```


<a id='Filtering-Records-1'></a>

# Filtering Records


```
Select rows within a given chromosome range. Provide chromosome range after this flag.
"--chromosome_range", "-r"
... -r chr1:20000-30000000 ...

Select rows within a given chromosome range. Provide chromosome range after this flag in format chr4:20000000-30000000
"--pass_filter", "-p"
... --pass_filter ...

Select variants matching list of chromosomal positions. Provide filename of text file formatted with two columns in .csv format.
"--positions_list", "-l"
... -l "example_positions_list.txt"
```


<a id='Selecting-and-Grouping-Samples-1'></a>

# Selecting and Grouping Samples


```
Select samples to include in visualization by providing tab delimited list of sample names (eg. samplenames.txt)
"--select_samples", "-x"
... -x example_list_of_sample_ids.txt

Group samples by common trait using user generated matrix key of traits and sample names following format guidelines in documentation. Provide file name of .csv file followed by trait to group by as it appears in the matrix key of traits and sample names.
"--group_samples", "-g"
... -g example_sample_traits_key.csv controls,cases ...
```

