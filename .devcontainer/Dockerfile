FROM braintwister/ubuntu-20.04:0.5

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
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
       software-properties-common \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN pip install \
      argparse \
      numpy \
      pandas \
      scipy \
      sklearn \
      emcee \
      lightkurve \
      loguru

RUN Rscript -e 'install.packages(c( \
      "argparser", \
      "devtools", \
      "numDeriv", \
      "lomb", \
      "tidyverse", \
      "modelr" \
    ), dependencies=TRUE)'

RUN Rscript -e 'devtools::install_version("ifultools", "2.0-23")'