FROM ubuntu:16.04

MAINTAINER Kamil Kwiek <kamil.kwiek@continuum.io>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    git mercurial subversion

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENV PATH /opt/conda/bin:$PATH

# ----------------
# modflow
# ----------------

RUN apt-get -y install \
    wget \
    unzip \
    git \
    gfortran \
    gcc \
    && apt-get clean \
    && rm -rf /tmp/* \
    && conda install numpy=1.11.2=py27_0 -y

RUN wget http://water.usgs.gov/ogw/modflow-nwt/MODFLOW-NWT_1.1.2.zip \
    && unzip MODFLOW-NWT_1.1.2.zip \
    && pip install https://github.com/modflowpy/pymake/zipball/master \
    && python -c "import os, pymake; pymake.main(srcdir='MODFLOW-NWT_1.1.2/src/', target='mfnwt', fc='gfortran', cc='gcc', makeclean=True, expedite=False, dryrun=False, double=False, debug=False)" \
    && rm -r MODFLOW-NWT_1.1.2 \
    && rm MODFLOW-NWT_1.1.2.zip

# ----------------
# HydroModelBuilder
# ----------------
RUN conda install conda-build
RUN mkdir /app
RUN git clone https://github.com/daniel-partington/HydroModelBuilder.git /app/HydroModelBuilder
RUN cd /app/HydroModelBuilder; git reset --hard a9221513feda93d35530ced0bcf1ff48fd570f1c
RUN conda build /app/HydroModelBuilder/docker_build/conda.recipe; cd /app

RUN git clone https://github.com/daniel-partington/CampaspeModel.git  /app/CampaspeModel
RUN cd /app/CampaspeModel; git reset --hard 6332960687118bab3a194fe62934928798639268
RUN conda build /app/CampaspeModel/docker_build/conda.recipe; cd /app

# change backend as "windows" don't play with docker. http://matplotlib.org/users/customizing.html#the-matplotlibrc-file
RUN MATPLOTLIB_RC_LOCATION=$(python -c "import matplotlib; print matplotlib.matplotlib_fname()") \
    && sed -i.bak 's/backend      : tkagg/backend      : Agg/gi' $MATPLOTLIB_RC_LOCATION
