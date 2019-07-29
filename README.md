# VIVA: A VCF File Visualization Tool and VariantVisualization.jl

## Visualization of Genomic Variants from VCF Files


| MacOS / Linux | Windows | License | Test Coverage | Documentation | Lifecycle |
| --- | ---- | ------ | ------ | ---- | ---- |
|[![Travis](https://img.shields.io/travis/compbiocore/VariantVisualization.jl/master.svg?style=flat-square)](https://travis-ci.org/compbiocore/VariantVisualization.jl)|[![Build status](https://ci.appveyor.com/api/projects/status/67hyn6rckulwr2dj/branch/master?svg=true)](https://ci.appveyor.com/project/fernandogelin/variantvisualization-jl/branch/master)|[![License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](https://github.com/compbiocore/VariantVisualization.jl/blob/master/LICENSE.md)|[![Coverage Status](https://coveralls.io/repos/github/compbiocore/VariantVisualization.jl/badge.svg?branch=master)](https://coveralls.io/github/compbiocore/VariantVisualization.jl?branch=master)|[![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/stable) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/latest) | ![Lifecycle](https://img.shields.io/badge/lifecycle-active-green.svg?style=flat-square) |


## Overview

VariantVisualization.jl is a package we built specifically to power the genetics visualization tool, *VIVA*.

*VIVA* is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files. It has been designed for clinicians and bioinformaticians to explore their VCF files visually. In a single command, users can extract genotype or read depth information and plot trends in interactive categorical heatmaps and scatter plots of average read depth values. VIVA offers a robust set of filters to select variants and samples of interest for analysis. VIVA is especially useful in early data exploration for identifying batch effect and sources of poor read depth in sequencing experiments, as well as identifying distribution of disease causing variants in a set of clinical samples.


## Getting Started:

## Installation

### Supported Operating Systems:

macOS ( Sierra, High Sierra, and Mojave ), Windows (7 and 10), and Linux.

### Step 1: Install Julia

Download [Julia]("https://julialang.org/downloads/") and install the language following the [platform specific instructions](https://julialang.org/downloads/platform.html).

Then, follow our [installation notes]("https://compbiocore.github.io/VariantVisualization.jl/latest/installation/") to add Julia to the path variable to run VIVA.

### Step 2: Install VariantVisualization.jl

To run the VIVA command line tool and VIVA Jupyter Notebook, you'll need to install our VariantVisualization.jl Julia package which powers VIVA.

To install VariantVisualization.jl:

1. Open the command line or PowerShell
2. Run the following block of code

>`julia`

>`]`

>`add VariantVisualization`

>`exit()`

### Step 3: Install the VIVA command line script

Download the VIVA tool script and save it to a working directory for your analysis. Save your VCF file in the working directory.

Copy and paste the following block of code into the command line or PowerShell:

>mkdir new_folder/

>cd new_folder/

>curl -L https://raw.githubusercontent.com/compbiocore/VariantVisualization.jl/master/viva > viva

### Optional Step: Install VIVA Jupyter Notebook

To install the VIVA Jupyer Notebook:

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).

Then, follow the in-notebook instructions to generate your plots.

## Run VIVA

Navigate in the Terminal or PowerShell to the directory containing the viva script run the VIVA command.

>cd new_folder/

>julia viva -f vcf.file arg1 arg2 arg3

We provide test files to run [EXAMPLES](https://compbiocore.github.io/VariantVisualization.jl/latest/examples/) after installation.

## For Developers

VIVA Jupyter notebook and the VIVA the command line tool are built with functions contained in our VariantVisualization.jl package.

Developers may contribute to these open source tools by using [functions contained within VariantVisualization.jl](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/) which are documented with docstrings.

We have included in-line comments within the code for the [VIVA command line tool](https://github.com/compbiocore/VariantVisualization.jl/tree/master/viva).

The VIVA Jupyter notebook is powered by a [main function](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/new_notebook_utils.jl) which takes arguments defined by the user in the notebook.


## Contributing and Questions

Contributions are welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/VariantVisualization.jl/issues
