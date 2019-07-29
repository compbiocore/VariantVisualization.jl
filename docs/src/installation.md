# Installation

### Supported Operating Systems:

macOS ( Sierra, High Sierra, and Mojave ), Windows (7 and 10), and Linux.

### Step 1: Install Julia

1. Download [Julia]("https://julialang.org/downloads/") and install the language following the [platform specific instructions](https://julialang.org/downloads/platform.html).

2. Then, follow add Julia to the path variable to run VIVA.

To add Julia to the PATH on Windows 7 or Windows 10:

Add the path to the Julia binaries (C:\Program Files\Julia\bin) to the PATH following the concise instructions [found here](https://www.java.com/en/download/help/path.xml)

To add Julia to the PATH on Mac run the following line in the Terminal:

> sudo ln -s /Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia

Be sure to replace "/Applications/Julia-1.1.app/..." to reflect the version of Julia you've downloaded.


*Linux Note*: To run on remote compute clusters, you may need to load the opengl and julia modules.


### Step 2: Install VariantVisualization.jl

To run the VIVA command line tool and VIVA Jupyter Notebook, you'll need to install our VariantVisualization.jl Julia package which powers VIVA.

To install VariantVisualization.jl:

1. Open the command line or PowerShell
2. Run the following block of code

```julia
julia
]
add VariantVisualization
exit()
```


### Step 3: Install the VIVA command line script

Download the VIVA tool script and save it to a working directory for your analysis. Save your VCF file in the working directory.

Copy and paste the following block of code into the command line or PowerShell:

>mkdir new_folder/

>cd new_folder/

>curl -L https://raw.githubusercontent.com/compbiocore/VariantVisualization.jl/master/viva > viva

### Optional Step: Install VIVA Jupyter Notebook

To install the VIVA Jupyer Notebook:

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).

Then, follow the in-notebook instructions to generate your plots.

### Latest Features

To stay up to date with cutting edge development features install VariantVisualization.jl from the Master branch.

from the Julia REPL (useful if using the PowerShell and don't have git installed):

```julia
]
add VariantVisualization#master
```

### For Developers

To add VariantVisualization in develop mode:

```julia
]
dev VariantVisualization
```

VIVA Jupyter notebook and the VIVA the command line tool are built with functions contained in our VariantVisualization.jl package.

Developers may contribute to these open source tools by using [functions contained within VariantVisualization.jl](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/) which are carefully documented with docstrings.

We have included in-line comments within the code for the [VIVA command line tool](https://github.com/compbiocore/VariantVisualization.jl/tree/master/viva).

The ***VIVA Jupyter notebook*** is powered by a [main function](https://github.com/compbiocore/VariantVisualization.jl/tree/master/src/new_notebook_utils.jl) which takes arguments defined by the user in the notebook. We welcome users to post in issues to request a new feature or bug fix.

## Installation Features Under Development

### Running VIVA with Docker or Docker Compose (Under Active Development)

Soon, you will be able to run VIVA using Docker images. This is not yet a supported feature. The instructions below will be helpful once this is supported.

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
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli:v0.3.9 /script/viva --save_remotely -f file.vcg -s pdf -o /data [...args]
```

- Example run:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli:v0.3.9 /script/viva --save_remotely -f file.vcf -s pdf -o /data [...args]
```

- On Windows:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli:v0.3.9 /script/viva --save_remotely -f file.vcf -s pdf -o /data [...args]
```

- Example run:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli:v0.3.9 /script/viva --save_remotely -f file.vcf -s pdf -o /data [...args]
```

##### Run the VIVA Jupyter Notebook from a Docker image:

Copy and run the following line from the terminal or Windows PowerShell:

- On Mac or Linux:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/notebook/data compbiocore/viva-notebook:v0.3.9
```

Go to the following url in your internet browser. You'll receive a token to enter into the url.

Go to `http://0.0.0.0:8888/?token=<enter token here>`

- On Windows:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "${pwd}":/home/jovyan/notebook/data compbiocore/viva-notebook:v0.3.9
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

-----
