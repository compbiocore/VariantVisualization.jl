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

------

## Using Docker and Docker Compose

If you don't want to install Julia and VariantVisualization, you can use the Docker images provided.
For that, first [install Docker](https://docs.docker.com/install/).

#### Using Docker

We provide two images, one with a Jupyter Notebook and one with a command line script for VIVA.

Create a project folder and navigate to it:
```shell
mkdir project_x
cd project_x
```

Make sure to add your project VCF files to that folder. That directory will be mapped to `/notebook/data` inside of the container.

Then, to run the Jupyter Notebook, from the terminal or Windows PowerShell:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/notebook/data compbiocore/viva-notebook
```
Go to `http://0.0.0.0:8888/?token=<enter token here>`

[Click here](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html) for more information about Jupyter Docker Images.

To run VIVA Command Line Tool:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli arg1 arg2 arg3
```

#### Using Docker Compose

To run the images with Docker Compose, copy the [`docker-compose.yml`](https://github.com/compbiocore/viva-docker/blob/master/docker-compose.yml) file to a local directory. From that same directory, run the commandas below.

!!! Note  
Your current directory will mount to `/notebook/data` in the notebook image and to `/data` in the CLI image.

- Notebook
```shell
docker-compose up viva-notebook
```

- Command Line Tool
```shell
docker-compose run viva -f /data/file.vcf arg2 arg3 ...
```

-----

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
