# ViVa.jl

#### Visualization of Variants


| MacOS / Linux | License | Test Coverage | Documentation | Lifecycle |
| --- | ---- | ------ | ------ | ---- |
|[![Travis](https://img.shields.io/travis/compbiocore/ViVa.jl/master.svg?style=flat-square)](https://travis-ci.org/compbiocore/ViVa.jl)| [![License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](https://github.com/compbiocore/ViVa.jl/blob/clean-up/LICENSE.md)| [![Codecov](https://img.shields.io/codecov/c/github/compbiocore/ViVa.jl.svg?style=flat-square)](https://codecov.io/gh/compbiocore/ViVa.jl/branch/master) | [![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/ViVa.jl/stable) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://compbiocore.github.io/ViVa.jl/latest) | ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg?style=flat-square) |

## Overview

ViVa.jl is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files and has been designed for clinicians and bioinformaticians to explore their VCF files visually. Users can quickly extract genotype or read depth information and plot trends in interactive categorical heatmaps and scatter plots of average read depth values. ViVa.jl offers a robust set of filters to select variants and samples of interest for analysis. ViVa.jl is especially useful in early data exploration for identifying batch effect and sources of poor read depth, as well as identifying distribution of disease causing variants in a set of clinical samples.


## Installation

Install ViVa.jl

```julia
using Pkg
Pkg.clone("https://github.com/compbiocore/ViVa.jl")
```

## Contributing and Questions

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/CAOS/issues
