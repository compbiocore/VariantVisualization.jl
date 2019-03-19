# ViVa.jl

#### Visualization of Variants


| MacOS / Linux | License | Test Coverage | Documentation | Lifecycle |
| --- | ---- | ------ | ------ | ---- |
|[![Travis](https://img.shields.io/travis/compbiocore/VIVA.jl/master.svg?style=flat-square)](https://travis-ci.org/compbiocore/VIVA.jl)| [![License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](https://github.com/compbiocore/VIVA.jl/blob/clean-up/LICENSE.md)| [![Codecov](https://img.shields.io/codecov/c/github/compbiocore/VIVA.jl.svg?style=flat-square)](https://codecov.io/gh/compbiocore/VIVA.jl/branch/master) | [![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/VIVA.jl/stable) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://compbiocore.github.io/VIVA.jl/latest) | ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg?style=flat-square) |

## Overview

VIVA.jl is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files and has been designed for clinicians and bioinformaticians to explore their VCF files visually. Users can quickly extract genotype or read depth information and plot trends in interactive categorical heatmaps and scatter plots of average read depth values. ViVa.jl offers a robust set of filters to select variants and samples of interest for analysis. ViVa.jl is especially useful in early data exploration for identifying batch effect and sources of poor read depth, as well as identifying distribution of disease causing variants in a set of clinical samples.


## Installation

### Command Line Tool

Add VIVA.jl in the Julia Pkg prompt.

### Jupyer Notebook

Install Jupyter and download the [VIVA Jupyter Notebook]().

### Latest Features

To stay up to date with cutting edge development features install VIVA.jl from the Master branch. 

Install ViVa.jl

```
git clone https://github.com/compbiocore/ViVa.jl
```

or

```julia
using Pkg
Pkg.clone("https://github.com/compbiocore/ViVa.jl")
```

## Contributing and Questions

Contributions are welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/VIVA.jl/issues
