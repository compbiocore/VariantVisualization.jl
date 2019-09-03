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

1. Download [Julia]("https://julialang.org/downloads/") and install the language following the [platform specific instructions](https://julialang.org/downloads/platform.html).

2. Then, follow add Julia to the path variable to run VIVA.

To add Julia to the PATH on Windows 7 or Windows 10:

Add the path to the Julia binaries (C:\Program Files\Julia\bin) to the PATH following the concise instructions [found here](https://www.java.com/en/download/help/path.xml)

To add Julia to the PATH on Mac run the following line in the Terminal:

> sudo ln -s /Applications/Julia-1.2.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia

Be sure to replace "/Applications/Julia-1.2.app/..." to reflect the version of Julia you've downloaded.


*Linux Note*: To run on remote compute clusters, you may need to load the opengl and julia modules.


### Step 2: Install VariantVisualization.jl

To run the VIVA command line tool and VIVA Jupyter Notebook, you'll need to install our VariantVisualization.jl Julia package which powers VIVA.

To install VariantVisualization.jl:

1. Open the command line or PowerShell
2. Run the following block of code

```julia
julia
]add VarianatVisualization
exit()
```
### Step 3: Run `viva`

#### Mac and Linux

On Mac and Linux, open another terminal window, navigate to your project folder and run:

```shell
viva -f filename.vcf -s <format> -o output/directory/
```

#### Windows

!!! Warning
    Viva will not work with Win32.

On windows, after installing VariantVisualization, open a new PowerShell and run:
```shell
viva -f filename.vcf -s <format> -o output/directory/
```

You'll then be prompted to select an application to open the script. Select the Julia executable, that is normally located
at `C:\Users\<username>\AppData\Local\Julia-<version>\bin\`.



### Optional Step: Install VIVA Jupyter Notebook

To install the VIVA Jupyer Notebook:

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).

Then, follow the in-notebook instructions to generate your plots.

### Latest Features

To stay up to date with cutting edge development features install VariantVisualization.jl from the Master branch.

From the Julia REPL:

```shell
julia
]add VariantVisualization#master
```

### For Developers

Install VariantVisualization in development mode:
```shell
julia
]dev VariantVisualization
```

VIVA Jupyter notebook and the VIVA the command line tool are built with functions contained in our VariantVisualization.jl package.

Developers may contribute to these open source tools by using [functions contained within VariantVisualization.jl](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/) which are carefully documented with docstrings.

We have included in-line comments within the code for the [VIVA command line tool](https://github.com/compbiocore/VariantVisualization.jl/tree/master/viva).

The ***VIVA Jupyter notebook*** is powered by a [main function](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/new_notebook_utils.jl) which takes arguments defined by the user in the notebook. We welcome users to post in issues to request a new feature or bug fix.

## Contributing and Questions

Contributions are welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/VariantVisualization.jl/issues
