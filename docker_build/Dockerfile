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

# change backend as "windows" don't play with docker. http://matplotlib.org/users/customizing.html#the-matplotlibrc-file
# RUN pip install \
#         matplotlib \
#         numpy \
#     && MATPLOTLIB_RC_LOCATION=$(python -c "import matplotlib; print matplotlib.matplotlib_fname()") \
#     && sed -i.bak 's/backend      : tkagg/backend      : Agg/gi' $MATPLOTLIB_RC_LOCATION \
#     && pip install flopy

# ----------------
# HydroModelBuilder
# ----------------
RUN mkdir app; cd app
RUN git clone https://github.com/daniel-partington/HydroModelBuilder.git
RUN cd HydroModelBuilder; git reset --hard 4da7b80873e0eedb25be11be849a26220cad402c; cd ..
RUN git clone https://github.com/daniel-partington/CampaspeModel.git
RUN cd CampaspeModel; git reset --hard 1e984b4d2d262a0912a19759e3170ef744de5698; cd ..
RUN conda install conda-build
RUN cd HydroModelBuilder; cd docker_build
RUN conda build conda.recipe
RUN cd ..; cd CampaspeModel; cd docker_build
RUN conda build conda.recipe


# ----------------
# CampaspeIM
# ----------------

# TODO:
#   Grab the code from respective bitbucket and github repos
#   run python install_packages.py for campaspeim
#   run python setup.py develop for each model
#   (the above should compile if necessary)

# docker run -it \
#     -v $PWD:/shared \
#     -e PYTHONPATH=/shared/CampaspeModel;/shared/HydroModelBuilder;/shared/CampaspeIM;/shared/CampasPyIHACRES \
#     campaspeim \
#     [run command]
