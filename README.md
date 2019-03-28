# VariantVisualization.jl

#### Visualization of Variants


| MacOS / Linux | License | Test Coverage | Documentation | Lifecycle |
| --- | ---- | ------ | ------ | ---- |
|[![Travis](https://img.shields.io/travis/compbiocore/VariantVisualization.jl/master.svg?style=flat-square)](https://travis-ci.org/compbiocore/VariantVisualization.jl)| [![License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](https://github.com/compbiocore/VariantVisualization.jl/blob/clean-up/LICENSE.md)| [![Codecov](https://img.shields.io/codecov/c/github/compbiocore/VariantVisualization.jl.svg?style=flat-square)](https://codecov.io/gh/compbiocore/VariantVisualization.jl/branch/master) | [![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/stable) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/latest) | ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg?style=flat-square) |

## Overview

VariantVisualization.jl is a package we built specifically to power the genetics visualization tool, *VIVA*.

*VIVA* is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files and has been designed for clinicians and bioinformaticians to explore their VCF files visually. Users can quickly extract genotype or read depth information and plot trends in interactive categorical heatmaps and scatter plots of average read depth values. VIVA offers a robust set of filters to select variants and samples of interest for analysis. VIVA is especially useful in early data exploration for identifying batch effect and sources of poor read depth, as well as identifying distribution of disease causing variants in a set of clinical samples.

To contribute to *VIVA*, developers may use the functions contained 


## Getting Started: *Installation*

### Supported Operating Systems:

#### macOS 
Mojave, High Sierra, and Mojave.

#### Windows
Windows 10, Windows 7

### Command Line Tool

1. Add VariantVisualization.jl in the Julia Pkg prompt.
2. Download the [VIVA](https://github.com/compbiocore/VariantVisualization.jl/tree/master/VIVA) tool script and save it to a working directory for your analysis.
3. Navigate to your working directory and follow the [VIVA manual](https://compbiocore.github.io/VariantVisualization.jl/latest/) to generate your plots.

### Jupyter Notebook

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/tree/master/VIVA.ipynb).
3. Follow the in-notebook instructions to generate your plots.

### Latest Features

To stay up to date with cutting edge development features install VariantVisualization.jl from the Master branch. 

Using git from the command line:

```
git clone https://github.com/compbiocore/VariantVisualization.jl
```

or from the Julia REPL (useful if using the PowerShell and don't have git installed):

```julia
using Pkg
Pkg.clone("https://github.com/compbiocore/VariantVisualization.jl")
```

### *For Developers*

VIVA Jupyter notebook and the VIVA the command line tool are built with functions contained in our VariantVisualization.jl package.

Developers may contribute to these open source tools by using [functions contained within VariantVisualization.jl](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/) which are carefully documented with docstrings.

We have included in-line comments within the code for the [VIVA command line tool](https://github.com/compbiocore/VariantVisualization.jl/tree/master/viva).

The ***VIVA Jupyter notebook*** is powered by a [main function](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/new_notebook_utils.jl) which takes arguments defined by the user in the notebook. We welcome users to post in issues to request a new feature or bug fix.


## Contributing and Questions

Contributions are welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/VariantVisualization.jl/issues
