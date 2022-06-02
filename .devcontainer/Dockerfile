FROM braintwister/ubuntu-20.04:0.5

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
 && add-apt-repository ppa:c2d4u.team/c2d4u4.0+ \
 && apt-get update \
 && export DEBIAN_FRONTEND=noninteractive \
 && apt-get -y install --no-install-recommends \
       apt-transport-https \
       build-essential \
       ca-certificates \
       dirmngr \
       gnupg \
       libblas-dev \
       libcurl4-gnutls-dev \
       liblapack-dev \
       libssl-dev \
       libxml2-dev \
       locales \
       python3-dev \
       python3-pip \
       r-base \
       r-base-dev \
       r-cran-devtools \
       r-cran-lomb \
       r-cran-modelr \
       r-cran-tidyverse \
       software-properties-common \
       tk \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c( \
      "argparser", \
      "numDeriv" \
    ), dependencies=TRUE)'

RUN Rscript -e 'devtools::install_version("ifultools", "2.0-23")'

RUN pip install \
      argparse \
      emcee \
      lightkurve \
      loguru \
      nbformat \
      numpy \
      pandas \
      plotly \
      scipy \
      sklearn

RUN pip install \
      black \
      pylint \
      pytest
