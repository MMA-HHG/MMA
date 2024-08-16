FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

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
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install -U --no-cache-dir \
    h5py==3.10.0 \
    jupyter==1.0.0 \
    jupyterlab==4.1.4 \
    matplotlib==3.8.3 \
    numpy==1.26.4
