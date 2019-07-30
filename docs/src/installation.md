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
]add VarianatVisualization
exit()
```
### Step 3: Run `viva`

#### Mac and Linux

On Mac and Linux, open another terminal window, navigate to your project folder and run:

```shell
viva -f filename.vcf -s <format> -o output/directory/
```

#### Windows

!!! Warning
    Viva will not work with Win32.

On windows, after installing VariantVisualization, open a new PowerShell and run:
```shell
viva -f filename.vcf -s <format> -o output/directory/
```

You'll then be prompted to select an application to open the script. Select the Julia executable, that is normally located
at `C:\Users\<username>\AppData\Local\Julia-<version>\bin\`.



### Optional Step: Install VIVA Jupyter Notebook

To install the VIVA Jupyer Notebook:

1. [Install Jupyter](https://jupyter.org/install)
2. Download the [VIVA Jupyter Notebook](https://github.com/compbiocore/VariantVisualization.jl/blob/master/VIVA.ipynb).

Then, follow the in-notebook instructions to generate your plots.

### Latest Features

To stay up to date with cutting edge development features install VariantVisualization.jl from the Master branch.

From the Julia REPL:

```shell
julia
]add VariantVisualization#master
```

### For Developers

Install VariantVisualization in development mode:
```shell
julia
]dev VariantVisualization
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

Make sure to add your project VCF files to that folder.


##### Run the VIVA Command Line Tool from a Docker image:

*Note*: Remember, you must use the flag `--save_remotely` when running VIVA by using Docker.

- On Mac or Linux:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli viva --save_remotely -f file.vcf -s pdf -o output
```

- Example run:
```shell
docker run -it --rm -v "$PWD":/data compbiocore/viva-cli viva --save_remotely -f file.vcf -s pdf -o output
```

- On Windows:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli viva --save_remotely -f file.vcf -s pdf -o output
```

- Example run:
```shell
docker run -it --rm -v "${pwd}":/data compbiocore/viva-cli viva --save_remotely -f file.vcf -s pdf -o output
```

##### Run the VIVA Jupyter Notebook from a Docker image:

Copy and run the following line from the terminal or Windows PowerShell:

- On Mac or Linux:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/notebook/data compbiocore/viva-notebook:v0.3.9
```

Go to the following url in your internet browser. You'll receive a token to enter into the url.

Go to `http://127.0.0.1:8888/?token=<enter token here>`

- On Windows:
```shell
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v "${pwd}":/home/jovyan/notebook/data compbiocore/viva-notebook:v0.3.9
```

Go to the following url in your internet browser. You'll receive a token to enter into the url.

Go to `http://127.0.0.1:8888/?token=<enter token here>`

[Click here](https://jupyter-docker-stacks.readthedocs.io/en/latest/index.html) for more information about Jupyter Docker Images.


-----
