
#Description


```
ViVa.jl is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files. ViVa.jl provides tools to quickly select variant records and samples to include in visualization.

ViVa.jl utilizes ArgParse.jl to parse command line arguments and GeneticVariation.jl to create a memory-efficient VCF Reader object which ViVa.jl utilizes to quickly extract variants of interest.

This documentation details each [REQUIRED] and optional argument, as well as necessary specific structure for user-generated input files for manipulating records and samples. Clearly named examples of all user-generated input files can be found in the /tests/ directory of the ViVa.jl repository.

ViVa.jl can be run through the command line with the following structure:

julia viva_cli.jl [--help] [COMMAND] [OPTIONS]

Example:
julia viva_cli.jl -v example.vcf -r chr1:20000-30000000 -s html --heatmap read_depth --line_chart samples
```


#Installation


> > > > > > > e6804e0d0fe65a23bcac9cfccb1395101d283ed2



Install ViVa.jl


To stay up to date with new features before release, please check out the master branch.


```julia
using Pkg
Pkg.clone("https://github.com/compbiocore/ViVa.jl")
```

