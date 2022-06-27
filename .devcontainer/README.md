Conda
=====

wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
bash ./Miniconda3-py38_4.12.0-Linux-x86_64.sh

conda env create --file conda-env.yml
conda env update --file conda-env.yml --prune
conda remove --name taco --all

conda activate taco
conda deactivate
