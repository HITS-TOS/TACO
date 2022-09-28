# Jupyterlab container

## Build docker container

Execute from TACO root directory:
```
docker build -t taco-jupyterlab -f .devcontainer/Dockerfile-jupyterlab .
```

## Run docker container

Insert data in <...> and execute from TACO root directory:

```
docker run -it --rm \
     -p <port>:8888 \
     --user root \
     -e NB_USER=<username> \
     -e NB_UID=<userid> \
     -e NB_GID=<groupid> \
     -e CHOWN_HOME=yes \
     -e CHOWN_HOME_OPTS="-R" \
     -w "/home/${NB_USER}" \
     -v "${PWD}":/home/<username>/work \
     taco-jupyterlab
```

Open the printed URL in your browser to access jupyterlab.

## Run unit tests (optional)

Open terminal in Jupyterlab and execute:

```
pip install pytest
export PYTHONPATH=./src
pytest
```


# Install TACO environment with conda

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
bash ./Miniconda3-py38_4.12.0-Linux-x86_64.sh

conda env create --file conda-env.yml
conda env update --file conda-env.yml --prune
conda remove --name taco --all

conda activate taco
conda deactivate
```


## Test pip packaging

# Build docker base container

```
docker build -t taco-base \
    --build-arg USERNAME=$(id -un) \
    --build-arg GROUPNAME=users \
    --build-arg USER_UID=$(id -u) \
    --build-arg USER_GID=$(id -g) \
    -f .devcontainer/Dockerfile-base .
```

# Run docker base container

```
docker run -it --rm \
    -v "${PWD}":/home/${USER}/work \
    taco-base bash
```

# Build python wheel

```
cd work
python3 -m build
```

# Test python wheel

```
pip install dist/taco-0.1-py3-none-any.whl
```
