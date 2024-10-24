[![Build Status](https://jenkins.h-its.org/buildStatus/icon?job=TOS%2FTACO%2Fmain)](https://jenkins.h-its.org/job/TOS/job/TACO/job/main/)

# Tools for Automated Characterisation of Oscillations (TACO)

The TACO modules will be restructured to be fully pythonic. Please find the heritage bash modules [here](README-legacy.md).


## Git usage

It is recommended to use git for downloading the TACO source code

```
git clone --recurse-submodules https://github.com/HITS-TOS/TACO.git
```

The dependency `sloscillations` is integrated as a git submodule and will be available using `--recurse-submodules` during git clone. If the flag was not used, it can be done afterwards with

```
git submodule update --init --recursive
```

Note: As long as the repository is private a [personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) is needed for the authentication.


## Tests (mostly not functioning, ignore)

Tests are implemented using pytest and can be executed with

```
python3 -m pytest
```


## Jupyterlab (functional, but not recommended)

The Jupyterlab docker container provides a comfortable way to perform TACO modules and can by started with

```
docker build -t taco-jupyterlab -f .devcontainer/Dockerfile-jupyterlab .
docker run -it --rm -p 8888:8888 taco-jupyterlab
```

Open the printed URL in your browser to access Jupyterlab. The jupyter notebook `work/pipeline.ipynb` is a good starting point.


## Install TACO with conda

Basing on [Miniconda](https://docs.conda.io/en/latest/miniconda.html), TACO can be installed with

```
conda env create
conda activate taco
```

## Install TACO with conda

Download and install the packages as per the requirements.txt file.

## Running high-throughput pipeline

For processing a long list of stars the high-throughput pipeline is available.
Before running the pipline, please execute
```
export PATH=$PWD/src:$PATH
export PYTHONPATH=$PWD/src:$PWD/libs/sloscillations:$PYTHONPATH
```
once from the TACO root directory.
Then the high-troughput pipline can be started with 
```
pipeline.py -i <input directory> -s <settings file>
```
taking every `<name>.dat` file in the `input directory` and write the results in a directory `<name>`.
A settings file with all entries is available at `pipeline/pipeline_settings_full.yaml`.

**Tip**

Copy the settings-file into a result directory and executing the pipline from there, leaves the run parameters documented.


## Tested operation system architectures

TACO docker-jupyterlab was tested on:
 - Linux (Ubuntu and CentOS)
 - MacOS (M1-Chip) (please consider https://github.com/HITS-TOS/TACO/issues/25)
 - Windows 11, Docker engine 4.15.0 using  WSL 2

TACO conda high-throughput pipeline was tested on:
 - Linux (CentOS)
