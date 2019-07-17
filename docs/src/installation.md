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
	* run `using Pkg`
	* run `Pkg.clone("https://github.com/compbiocore/VariantVisualization.jl")`
	* run `Pkg.instantiate()`
2. Download the [VIVA](https://github.com/compbiocore/VariantVisualization.jl/blob/master/viva) tool script and save it to a working directory for your analysis.
3. Navigate to your working directory and follow the [VIVA manual](https://compbiocore.github.io/VariantVisualization.jl/latest/) to generate your plots.

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

To run the images with Docker Compose, install Docker following the steps above and then install [Docker Compose[(https://docs.docker.com/compose/). Then copy the [docker-compose.yml](https://github.com/compbiocore/viva-docker/blob/master/docker-compose.yml) file to a local directory. From that same directory, run the command as it appears below.

*Note*: Your current directory will mount to `/notebook/data` in the notebook image and to `/data` in the CLI image.

- Notebook
```shell
docker-compose up viva-notebook
```

- Command Line Tool
```shell
docker-compose run viva -f --save_remotely /data/file.vcf arg2 arg3 ...
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
