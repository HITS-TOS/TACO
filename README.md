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


## Jupyterlab

The Jupyterlab docker container provides a comfortable way to perform TACO modules and can by started with

```
docker build -t taco-jupyterlab -f .devcontainer/Dockerfile-jupyterlab .
docker run -it --rm -p 8888:8888 taco-jupyterlab
```

Open the printed URL in your browser to access Jupyterlab. The jupyter notebook `work/pipeline.ipynb` is a good starting point.
