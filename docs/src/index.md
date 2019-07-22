![VIVA Logo](assets/VIVA_logo.png)

# Getting Started

# *VIVA Command Line Tool and Jupyter Notebook*

## Description

VIVA is a user-friendly command line tool built with our VariantVisualization.jl package for exploratory analysis and generation of publication quality graphics for variant analysis projects using Variant Call Format (VCF) files.

Variant selection and plotting is all executed in a single command.

We describe each of VIVA's arguments in this documentation under the Manual page.

VIVA is available as a Jupyter Notebook utility [here](https://github.com/compbiocore/VariantVisualization.jl/tree/master). Instructions for installing Jupyter, downloading VIVA Jupyter Notebook, and using the notebook are detailed in the *Jupyter Notebook* section of this documentation.

Formatting requirements for VIVA's input files are described in the Manual and clearly named examples of all user-generated input files can be found in the `/test/test_files` directory of the `VariantVisualization.jl` repository.

## General Use

To use VIVA, we recommend creating a new directory for storing your VCF file to analyze where output files will be saved. Alternatively, users may also provide paths to the VCF file and to preferred output file locations as command line arguments.

### Command Line
VIVA's general command line argument structure is as follows:

```
    julia viva -f file.vcf [OPTIONS]
```

From the command line or powershell, run the VIVA command line tool script which takes arguments from the command line and parses them with ArgParse.jl.

Example:

```
    julia viva -f example.vcf -r chr1:20000-30000000 -s pdf -m genotype,read_depth --avg_dp samples
```

To display a complete set of help instructions while using the tool, run VIVA with the help flag (`--help`, `-h`).

```
    julia viva -h
```

### Default Options:

By running VIVA with only a VCF filename:

```
    julia viva -f file.vcf
```

Default options will be used:

`--heatmap` = `genotype,read_depth`
`--save_format` = `html`
`--output_directory` = `output`
`--heatmap_title` = `vcf_filename`
`--y_axis_labels` = `chromosomes`
`--x_axis_labels` = `true`

These default settings generate a heatmap plots of genotype and read depth values of all variants for all sample ids within a VCF file.

We recommend using variant filters with most VCF files as there is too much data to plot or evaluate visually.

Specifically, we recommend visualizing fewer than 2000 variants at a time for effective visualization. However, VIVA uses memory efficient filtering and plotting and is capable of plotting >200,000 datapoints.

### Jupyter Notebook

Use the following steps to use the VIVA Jupyter Notebook utility:

1. Install Jupyter Notebook following the [platform specific instructions](https://plot.ly/python/ipython-notebook-tutorial/)

2. Download the [VIVA Jupyter Notebook](PATH) to a working directory containing your VCF file.

3. Open the Julia REPL on the command line from any directory.

4. Run `using IJulia` and then `notebook()`

5. Navigate to the directory containing the VIVA Jupyter Notebook *VIVA.ipynb* and double click to open.

6. Follow the step-by-step instructions within the notebook to generate your figures.

### Running VIVA with Docker or Docker Compose

Alternatively, you can run VIVA using the Docker images we've provided if you don't want to install Julia and the VariantVisualization.jl Julia package. You may only save images to HTML format using the Docker, for now, due to technical limitations of dependency packages. We've actively developing a feature to save to all formats using Docker.

To run VIVA from a Docker image, first [install Docker](https://docs.docker.com/install/).

Then double-click the Docker.app in the Applications folder to start Docker. You will see a whale icon in the top status bar to indicate that Docker is running and accessible from the terminal. You can quit Docker once you are finished using VIVA by clicking the Docker whale icon in the top status bar and clicking "Quit Docker Desktop."

#### Using Docker

*Note*: You must use the flag `--save_remotely` when running VIVA by using Docker.

Once Docker is running, you can run VIVA by running the Docker commands below in the Mac/Linux terminal or Windows PowerShell.

We provide two images, one with a Jupyter Notebook and one with a command line script for VIVA. You can run VIVA in a single command using these images. The command consists of calls to run the Docker image followed by the usual VIVA options.

To run the images, follow these steps:

Create a project folder and navigate to it:
```shell
mkdir project_x
cd project_x
```

Make sure to add your project VCF files to that folder. That directory will be mapped to `/notebook/data` inside of the container.

When entering the filename of the VCF file and files to support filtering options, you should include `/data/...` in the path to your files.

##### Run the VIVA Command Line Tool from a Docker image:

*Note*: Remember, you must use the flag `--save_remotely` when running VIVA by using Docker.

- On Mac or Linux:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli:v0.3.8 --save_remotely arg1 arg2 arg3
```

- Example run:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli:v0.3.8 --save_remotely -f file.vcf -p
```

- On Windows:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli:v0.3.8 --save_remotely arg1 arg2 arg3
```

- Example run:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli:v0.3.8 --save_remotely -f file.vcf -p
```

##### Run the VIVA Jupyter Notebook from a Docker image:

Copy and run the following line from the terminal or Windows PowerShell:

- On Mac or Linux:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/notebook/data compbiocore/viva-notebook
```

Go to the following url in your internet browser. You'll receive a token to enter into the url.

Go to `http://0.0.0.0:8888/?token=<enter token here>`

- On Windows:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "${pwd}":/home/jovyan/notebook/data compbiocore/viva-notebook
```

Go to the following url in your internet browser. You'll receive a token to enter into the url.

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
docker-compose run viva -f file.vcf --save_remotely arg3 arg4 ...
```

## Continue reading for:

* [Variant and Sample Selection](https://compbiocore.github.io/VariantVisualization.jl/stable/filtering_vcf/)

* [Plotting Options](https://compbiocore.github.io/VariantVisualization.jl/stable/plotting/)
