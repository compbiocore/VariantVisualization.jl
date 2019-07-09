# VariantVisualization.jl
## Visualization of Variants


| MacOS / Linux | Windows | License | Test Coverage | Documentation | Lifecycle |
| --- | ---- | ------ | ------ | ---- | ---- |
|[![Travis](https://img.shields.io/travis/compbiocore/VariantVisualization.jl/master.svg?style=flat-square)](https://travis-ci.org/compbiocore/VariantVisualization.jl)|[![Build status](https://ci.appveyor.com/api/projects/status/67hyn6rckulwr2dj/branch/master?svg=true)](https://ci.appveyor.com/project/fernandogelin/variantvisualization-jl/branch/master)|[![License](https://img.shields.io/badge/license-MIT-orange.svg?style=flat-square)](https://github.com/compbiocore/VariantVisualization.jl/blob/master/LICENSE.md)|[![Coverage Status](https://coveralls.io/repos/github/compbiocore/VariantVisualization.jl/badge.svg?branch=master)](https://coveralls.io/github/compbiocore/VariantVisualization.jl?branch=master)|[![Docs](https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/stable) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square)](https://compbiocore.github.io/VariantVisualization.jl/latest) | ![Lifecycle](https://img.shields.io/badge/lifecycle-active-green.svg?style=flat-square) |


## Overview

VariantVisualization.jl is a package we built specifically to power the genetics visualization tool, *VIVA*.

*VIVA* is a user-friendly command line tool for creating publication quality graphics from Variant Call Format (VCF) files. It has been designed for clinicians and bioinformaticians to explore their VCF files visually. In a single command, users can extract genotype or read depth information and plot trends in interactive categorical heatmaps and scatter plots of average read depth values. VIVA offers a robust set of filters to select variants and samples of interest for analysis. VIVA is especially useful in early data exploration for identifying batch effect and sources of poor read depth, as well as identifying distribution of disease causing variants in a set of clinical samples.


## Getting Started:

Note: Once you have set up VIVA, you can quickly run the command line tool [EXAMPLES](https://compbiocore.github.io/VariantVisualization.jl/latest/examples/) found in the documentation.


## Installation

### Supported Operating Systems:

macOS ( Sierra, High Sierra, and Mojave ), Windows, and Linux.

To use VIVA, you must download the Julia programming language version >=1.0 and install the VariantVisualization.jl Julia package as well as the VIVA script.

Expected Time for Installation: Installation time depends on your network bandwidth, but should take less than 10 minutes for VIVA installation to install all dependency packages. Installing and using Julia packages for the first time takes longer than when using them in subsequent sessions.

*Note*: When installing VariantVisualization.jl and running VIVA remote compute clusters, you may need to load the OpenGl module in addition to loading the Julia module.

### Command Line Tool

1. Add VariantVisualization.jl using Pkg in the Julia REPL:
	a. Open the Julia REPL by typing `julia` into the command line
	b. Enter the Pkg manager by entering ']' into the REPL
	c. Enter `add VariantVisualization` in the Pkg manager. This will install all of VIVA's dependencies.
2. Download the [VIVA](https://github.com/compbiocore/VariantVisualization.jl/blob/master/viva) tool script and save it to a working directory for your analysis.
3. Navigate to your working directory and follow the [VIVA manual](https://compbiocore.github.io/VariantVisualization.jl/stable/) to generate your plots.

### Jupyter Notebook

1. [Install Jupyter](https://jupyter.org/install)
2. Install the VariantVisualization.jl Julia package following the Command Line Tool installation instructions above.
3. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).
4. Follow the in-notebook instructions to generate your plots.

### Running VIVA with Docker or Docker Compose

Alternatively, you can run VIVA using the Docker images we've provided if you don't want to install Julia and the VariantVisualization.jl Julia package.

To run VIVA from a Docker image, first [install Docker](https://docs.docker.com/install/).

Then double-click the Docker.app in the Applications folder to start Docker. You will see a whale icon in the top status bar to indicate that Docker is running and accessible from the terminal. You can quit Docker once you are finished using VIVA by clicking the Docker whale icon in the top status bar and clicking "Quit Docker Desktop."

#### Using Docker

*Note* You must use the flag `--save_remotely` when running VIVA by using Docker.

Once Docker is running, you can run VIVA by running the Docker commands below in the Mac/Linux terminal or Windows PowerShell.

We provide two images, one with a Jupyter Notebook and one with a command line script for VIVA. You can run VIVA in a single command using these images. The command consists of calls to run the Docker image followed by the usual VIVA options.

To run the images, follow these steps:

Create a project folder and navigate to it:
```shell
mkdir project_x
cd project_x
```

Make sure to add your project VCF files to that folder. That directory will be mapped to `/notebook/data` inside of the container.

##### Run the VIVA Command Line Tool from a Docker image:

*Note* Remember, you must use the flag `--save_remotely` when running VIVA by using Docker.

- On Mac or Linux:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli --save_remotely arg1 arg2 arg3
```

- Example run:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli --save_remotely -f file.vcf -p -s pdf
```

- On Windows:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli --save_remotely arg1 arg2 arg3
```

- Example run:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli --save_remotely -f file.vcf -p -s pdf
```

##### Run the VIVA Jupyter Notebook from a Docker image:

Copy and run the following line from the terminal or Windows PowerShell:

- On Mac or Linux:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/notebook/data compbiocore/viva-notebook
```
Go to `http://0.0.0.0:8888/?token=<enter token here>`

- On Windows:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "${pwd}":/notebook/data compbiocore/viva-notebook
```
Go to `http://0.0.0.0:8888/?token=<enter token here>`

[Click here](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html) for more information about Jupyter Docker Images.

#### Using Docker Compose

To run the images with Docker Compose, copy the [`docker-compose.yml`](https://github.com/compbiocore/viva-docker/blob/master/docker-compose.yml) file to a local directory. From that same directory, run the command as it appears below.

*Note*: Your current directory will mount to `/notebook/data` in the notebook image and to `/data` in the CLI image.

- Notebook
```shell
docker-compose up viva-notebook
```

- Command Line Tool
```shell
docker-compose run viva -f --save_remotely /data/file.vcf arg2 arg3 ...
```

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


### For Developers

VIVA Jupyter notebook and the VIVA the command line tool are built with functions contained in our VariantVisualization.jl package.

Developers may contribute to these open source tools by using [functions contained within VariantVisualization.jl](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/) which are carefully documented with docstrings.

We have included in-line comments within the code for the [VIVA command line tool](https://github.com/compbiocore/VariantVisualization.jl/tree/master/viva).

The ***VIVA Jupyter notebook*** is powered by a [main function](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/new_notebook_utils.jl) which takes arguments defined by the user in the notebook. We welcome users to post in issues to request a new feature or bug fix.


## Contributing and Questions

Contributions are welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[issues-url]: https://github.com/compbiocore/VariantVisualization.jl/issues
