#Installation

##Install Julia v1.1.0
Download [Julia]("https://julialang.org/downloads/")

##Install the VariantVisualization.jl
Add `VariantVisualization.jl` at the package prompt in the Julia v1.1 REPL.

```
julia

*press the ']' key to enter the package prompt*

(v1.1) pkg> add VariantVisualization

```

If successful, `VariantVisualization.jl` will be installed in your Julia packages directory (.julia/packages/) and VIVA command line tool should install with this. Adding `VariantVisualization.jl` runs a script to create an alias for the command line tool script in the bash shell. This allows the user to call VIVA from any directory.


##Install the Jupyter Notebook

If you plan to use the Jupyter Notebook VIVA utility, install Jupyter then download the [VIVA Notebook]().

If you already have Jupyter installed, update following [these instructions](https://jupyter.readthedocs.io/en/latest/projects/upgrade-notebook.html)


##New Features
To stay up to date with new features before official version release, please check out the master branch.

```
julia

using Pkg
Pkg.clone("https://github.com/compbiocore/VariantVisualization.jl")
```
