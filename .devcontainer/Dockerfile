FROM ubuntu:20.04

LABEL maintainer="Bernd Doser <bernd.doser@h-its.org>"

ARG USERNAME=user
ARG GROUPNAME=users
ARG USER_UID=1000
ARG USER_GID=1000
ARG TIMEZONE='Europe/Berlin'
ARG DEBIAN_FRONTEND=noninteractive

RUN echo $TIMEZONE > /etc/timezone && \
  apt-get update && apt-get install -y tzdata && \
  rm /etc/localtime && \
  ln -snf /usr/share/zoneinfo/$TIMEZONE /etc/localtime && \
  dpkg-reconfigure -f noninteractive tzdata && \
  apt-get clean

RUN apt-get update \
 && apt-get -y install \
       apt-transport-https \
       build-essential \
       ca-certificates \
       curl \
       dirmngr \
       git \
       gnupg \
       gpg-agent \
       libblas-dev \
       libcurl4-gnutls-dev \
       libfontconfig1-dev \
       libfreetype6-dev \
       libfribidi-dev \
       libharfbuzz-dev \
       libjpeg-dev \
       liblapack-dev \
       libpng-dev \
       libssl-dev \
       libtiff5-dev \
       libxml2-dev \
       locales \
       make \
       ninja-build \
       python3 \
       python3-dev \
       python3-pip \
       python3-setuptools \
       python3-tk \
       r-base \
       r-base-core \
       r-base-dev \
       software-properties-common \
       sudo \
       tk \
       vim \
       wget \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c( \
      "argparser", \
      "languageserver", \
      "numDeriv", \
      "devtools", \
      "lomb", \
      "modelr", \
      "tidyverse" \
    ), dependencies=TRUE)'

RUN Rscript -e 'devtools::install_version("ifultools", "2.0-23")'

RUN /usr/bin/python3 -m pip install --upgrade pip

RUN pip install \
      argparse \
      emcee \
      lightkurve \
      loguru \
      nbformat \
      numba \
      numpy \
      pandas \
      plotly \
      pyyaml \
      rpy2==3.3.* \
      scipy \
      sklearn

RUN pip install \
      black \
      ipykernel \
      pylint \
      pytest \
      pytest-timeout

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

RUN useradd --uid ${USER_UID} -m ${USERNAME} \
 && adduser ${USERNAME} sudo \
 && echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

USER ${USERNAME}
