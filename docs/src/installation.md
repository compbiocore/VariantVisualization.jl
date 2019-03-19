#Installation

##Install Julia v1.1
Follow the [instructions]("") for installing Julia v1.1 on your system.

##Install the VIVA.jl
Add `VIVA.jl` at the package prompt in the Julia v1.1 REPL.

```
julia

*press the ']' key to enter the package prompt*

(v1.1) pkg> add VIVA

```

If successful, `VIVA.jl` will be installed in your Julia packages directory (.julia/packages/) and VIVA command line tool should install with this. Adding `VIVA.jl` runs a script to create an alias for the command line tool script in the bash shell. This allows the user to call VIVA from any directory.


##Install the Jupyter Notebook

If you plan to use the Jupyter Notebook VIVA utility, install Jupyter then download the [VIVA Notebook]().

***Fernando, let's decide how to make VIVA notebook accessible***

If you already have Jupyter installed, update following [these instructions](https://jupyter.readthedocs.io/en/latest/projects/upgrade-notebook.html)


##New Features
To stay up to date with new features before official version release, please check out the master branch.

***Fernando, rather than using Pkg, do we recommend using the packages prompt ']' ?***

```
julia

using Pkg
Pkg.clone("https://github.com/compbiocore/VIVA.jl")
```
