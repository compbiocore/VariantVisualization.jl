#Installation

### Install Julia v1.1.0
Download [Julia]("https://julialang.org/downloads/")

### Supported Operating Systems:

#### macOS 

Sierra, High Sierra, and Mojave.

#### Windows

Windows 10, Windows 7.

#### Linux

*Note*: To run on remote compute clusters, you may need to load opengl module along with julia/1.1.0.

### Command Line Tool

1. Add VariantVisualization.jl using Pkg in the Julia REPL:
	run `using Pkg`
	run `Pkg.clone("https://github.com/compbiocore/VariantVisualization.jl")`
	run `Pkg.instantiate()`
2. Download the [VIVA](https://github.com/compbiocore/VariantVisualization.jl/blob/master/viva) tool script and save it to a working directory for your analysis.
3. Navigate to your working directory and follow the [VIVA manual](https://compbiocore.github.io/VariantVisualization.jl/latest/) to generate your plots.

### Jupyter Notebook

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).
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

