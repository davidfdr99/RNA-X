FROM --platform=linux/amd64 ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# File Author / Maintainer
MAINTAINER David Fandrei <david.fandrei@gustaveroussy.fr>

USER root

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update \
 && apt-get install -y \
  apt-utils \
  sudo \
  git \
  less \
  wget \
  tree \
  graphviz \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
  apt-get install -y openjdk-8-jdk && \
  apt-get install -y ant && \
  apt-get clean;

RUN apt-get install ca-certificates-java && \
  update-ca-certificates -f;

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

ENV PATH="/opt/conda/bin:$PATH"

RUN conda install -c conda-forge mamba

COPY ../environment.yml .

RUN /bin/bash -c "time mamba env create -f environment.yml"

EXPOSE 5003

RUN echo -e "#! /bin/bash\n\n# script to activate the conda environment" > ~/.bashrc \
  && echo "export PS1='Docker> '" >> ~/.bashrc \
  && conda init bash \
  && echo "\nconda activate rnax" >> ~/.bashrc \
  && conda clean -a

# environment variables
ENV LANG en_US.utf-8
ENV BASH_ENV ~/.bashrc


