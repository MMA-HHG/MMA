FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

ENV MSM_PATH=/MMA
ENV MULTISCALE_WORK_DIR=$MSM_PATH/work_dir
ENV MULTISCALE_DEMOS=$MULTISCALE_WORK_DIR

ENV PYTHONPATH=$PYTHONPATH:$MSM_PATH/shared_python

ENV CUPRAD_HOME=$MSM_PATH/CUPRAD
ENV CUPRAD_BUILD=$MSM_PATH/CUPRAD/build
ENV CUPRAD_SCRIPTS=$MSM_PATH/CUPRAD/scripts
ENV CUPRAD_PYTHON=$MSM_PATH/CUPRAD/python
ENV PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON

ENV TDSE_1D_HOME=$MSM_PATH/1DTDSE
ENV TDSE_1D_PYTHON=$MSM_PATH/1DTDSE/python
ENV TDSE_1D_SCRIPTS=$MSM_PATH/1DTDSE/scripts
ENV TDSE_1D_SLURM=$MSM_PATH/1DTDSE/slurm
ENV TDSE_1D_BUILD=$MSM_PATH/1DTDSE/build
ENV PYTHONPATH=$PYTHONPATH:$TDSE_1D_HOME

ENV HANKEL_HOME=$MSM_PATH/Hankel

ENV MULTISCALE_HOME=$MSM_PATH
ENV MULTISCALE_SCRIPTS=$MSM_PATH/multiscale/scripts

ENV JUPYTER_EXAMPLES=$MSM_PATH/jupyter_examples

ENV FSPA_PATH=$MSM_PATH/FSPA

ENV CC=mpicc
ENV FC=mpifort

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        gfortran \
        make \
        cmake \
        mpich \
        libfftw3-dev \
        libfftw3-mpi-dev \
        libhdf5-openmpi-dev \
        python3 \
        python3-pip \
        python3-setuptools \
        python3-wheel \
        git \
        ffmpeg \
        bash-completion \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir \
    h5py==3.10.0 \
    jupyter==1.0.0 \
    jupyterlab==4.1.4 \
    matplotlib==3.8.3 \
    numpy==1.26.4 \
    scipy==1.11.3 \ 
    ipykernel==6.29.5 \
    multiprocess==0.70.17 \
    jupyterlab-h5web==12.5.0 \
    git+https://github.com/vabekjan/universal_input.git@pipable#egg=create_universal_HDF5 

# 3. Make JupyterLab terminals start Bash (login shell) instead of sh
RUN jupyter server --generate-config && \
    printf "\nc.ServerApp.terminado_settings = {'shell_command': ['/bin/bash','-l']}\n" \
    >> /root/.jupyter/jupyter_server_config.py

ENV SHELL=/bin/bash