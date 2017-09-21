FROM ubuntu:16.04

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

ENTRYPOINT [ "/usr/bin/tini", "--" ]

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
    && conda install numpy=1.11.2=py27_0 matplotlib gdal xlrd pysal simplejson pytables -c conda-forge -y

# Have to manually set GDAL_DATA env variable
ENV GDAL_DATA /opt/conda/share/gdal

RUN pip install pyevtk more_itertools jsmin

RUN wget http://water.usgs.gov/ogw/modflow-nwt/MODFLOW-NWT_1.1.2.zip \
    && unzip MODFLOW-NWT_1.1.2.zip \
    && pip install https://github.com/modflowpy/pymake/zipball/master \
    && python -c "import os, pymake; pymake.main(srcdir='MODFLOW-NWT_1.1.2/src/', target='mfnwt', fc='gfortran', cc='gcc', makeclean=True, expedite=False, dryrun=False, double=False, debug=False)" \
    && rm -r MODFLOW-NWT_1.1.2 \
    && rm MODFLOW-NWT_1.1.2.zip


# change backend as "windows" don't play with docker. http://matplotlib.org/users/customizing.html#the-matplotlibrc-file
RUN MATPLOTLIB_RC_LOCATION=$(python -c "import matplotlib; print matplotlib.matplotlib_fname()") \
    && sed -i.bak 's/backend      : qt5agg/backend      : Agg/gi' $MATPLOTLIB_RC_LOCATION


# -----------------
# HydroModelBuilder
# -----------------
ENV CONDA_PY 2.7
RUN pip install git+https://github.com/daniel-partington/HydroModelBuilder@32f8bb6bc8885b3ac4685f4459d772220c242b30

# while in development install this at runtime ie. PYTHONPATH to shared folder

# ----------------
# CampaspeModel
# ----------------
# RUN mkdir /app. \
#     && git clone https://github.com/daniel-partington/CampaspeModel.git /app/CampaspeModel \
#     && cd /app/CampaspeModel; git reset --hard 0a2829992aaab4615bf3d147a5d038cf6cbc276d

# Add CampaspeModel to PYTHONPATH
# ENV PYTHONPATH /app/CampaspeModel:$PYTHONPATH


# To use CampaspeModel, you must have a shared directory with model data set up, and
# the relevant settings in model_config.json must point to this directory.
# This defaults to /shared in the model_config.json file.

# Docker run command
# docker run -v [path to share folder on host]:[docker directory] -it [container]
# so, for example
# docker run -v C:/shared_folder:/shared -it campaspemodel
