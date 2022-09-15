[![Build Status](https://jenkins.h-its.org/buildStatus/icon?job=TOS%2FTACO%2Fmain)](https://jenkins.h-its.org/job/TOS/job/TACO/job/main/)

# Tools for Automated Characterisation of Oscillations (TACO)

The TACO modules will be restructured to be fully pythonic. Please find the heritage bash modules [here](README-legacy.md).


## Jupyterlab

The Jupyterlab docker container provides a comfortable way to perform TACO modules and can by started with

```
docker build -t taco-jupyterlab -f .devcontainer/Dockerfile-jupyterlab .
docker run -it --rm \
     -p 8888:8888 \
     --user root \
     -e NB_USER=$(id -un) \
     -e NB_UID=$(id -u) \
     -e NB_GID=$(id -g) \
     -e CHOWN_HOME=yes \
     -e CHOWN_HOME_OPTS="-R" \
     -w "/home/${NB_USER}" \
     -v "${PWD}":/home/$(id -un)/work \
     taco-jupyterlab
```

Open the printed URL in your browser to access Jupyterlab. The jupyter notebook `juputer/pipeline.ipynb` is a good starting point.
