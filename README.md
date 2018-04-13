# ViVa

[![Build Status](https://travis-ci.org/compbiocore/ViVa.jl.svg?branch=master)](https://travis-ci.org/compbiocore/ViVa.jl)

[![Coverage Status](https://coveralls.io/repos/compbiocore/ViVa.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/compbiocore/ViVa.jl?branch=master)

[![codecov.io](http://codecov.io/github/compbiocore/ViVa.jl/coverage.svg?branch=master)](http://codecov.io/github/compbiocore/ViVa.jl?branch=master)


ViVa.jl is a package for parsing and visualizing field information contained within VCF files. Users may apply filters to variant rows and sample columns of the VCF file before visualizing data from a chosen vcf data field.

## Command line usage

### ARGUMENTS Guide

ARGS[1] = VCF file
ARGS[2] = format to save graphic as (pdf, svg, etc.)
ARGS[3] = field to visualize (-gt, -dp)
ARGS[4] = variant filter to apply (-l, -r, -a)
ARGS[5] = PASS FILTER only
ARGS[6] = reorder_columns
ARGS[7] = phenotype matrix filename for reordering columns (in .csv)
ARGS[8] = chromosome range, significant variant list filename, or nothing (.)
ARGS[9] = select_columns
ARGS[10] = filename of list of sample names to select
ARGS[11] = phenotype to sort columns by (used by ARGS[7]) (e.g. case_control_status)

#### Example 1:  

This is a command line example for visualizing genotype of variants which passed QC filters within a chromosome range across a selected group of samples grouped by case_control_status:

```julia
julia masterv0.1.jl file.vcf pdf -gt -r pass_only reorder_columns phenotype_matrix.csv chr1:10000000-15000000 select_columns sample_names.tsv case_control_status
```

### Example 2: 
This is a command line example for visualizing read depth for all variants with no filters applied:

```julia
julia masterv0.1.jl file.vcf pdf -d -a . . . . . . .
```

### Example 2: 
This is a command line example for visualizing read depth for all variants with sample columns grouped by case_control_status

```julia
julia masterv0.1.jl file.vcf pdf -d -a . reorder_columns phenotype_matrix.csv . . . case_control_status
```