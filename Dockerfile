FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

# Default container user. Linux users can override UID/GID at build time.
ARG USERNAME=active_electron
ARG USER_UID=1000
ARG USER_GID=1000

# Project paths used by MMA scripts and examples.
ENV MMA_PATH=/MMA
ENV MULTISCALE_WORK_DIR=/MMA/work_dir
ENV MULTISCALE_DEMOS=/MMA/work_dir

ENV CUPRAD_HOME=/MMA/CUPRAD
ENV CUPRAD_BUILD=/MMA/CUPRAD/build
ENV CUPRAD_SCRIPTS=/MMA/CUPRAD/scripts
ENV CUPRAD_PYTHON=/MMA/CUPRAD/python

ENV TDSE_1D_HOME=/MMA/1DTDSE
ENV TDSE_1D_PYTHON=/MMA/1DTDSE/python
ENV TDSE_1D_SCRIPTS=/MMA/1DTDSE/scripts
ENV TDSE_1D_SLURM=/MMA/1DTDSE/slurm
ENV TDSE_1D_BUILD=/MMA/1DTDSE/build

ENV HANKEL_HOME=/MMA/Hankel
ENV MULTISCALE_HOME=/MMA
ENV MULTISCALE_SCRIPTS=/MMA/multiscale/scripts
ENV JUPYTER_EXAMPLES=/MMA/jupyter_examples
ENV FSPA_PATH=/MMA/FSPA

ENV PYTHONPATH=/MMA/shared_python:/MMA/CUPRAD/python:/MMA/1DTDSE

# Use MPI compiler wrappers by default.
ENV CC=mpicc
ENV FC=mpifort
ENV SHELL=/bin/bash

# System dependencies for native code, Python, Jupyter, and development.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
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
        sudo \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies are kept in the standard pip requirements format.
COPY docker/requirements.txt /tmp/requirements.txt
RUN python3 -m pip install --upgrade --no-cache-dir pip \
    && python3 -m pip install --no-cache-dir -r /tmp/requirements.txt \
    && rm -f /tmp/requirements.txt

# Make JupyterLab terminals start Bash instead of sh.
RUN mkdir -p /etc/jupyter \
    && printf "c.ServerApp.terminado_settings = {'shell_command': ['/bin/bash', '-l']}\n" \
        > /etc/jupyter/jupyter_server_config.py

# Create a normal user for MPI runs; sudo keeps the image useful for development.
RUN groupadd --gid "${USER_GID}" "${USERNAME}" \
    && useradd \
        --uid "${USER_UID}" \
        --gid "${USER_GID}" \
        --create-home \
        --shell /bin/bash \
        "${USERNAME}" \
    && install -d -m 0755 /etc/sudoers.d \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > "/etc/sudoers.d/${USERNAME}" \
    && chmod 0440 "/etc/sudoers.d/${USERNAME}" \
    && mkdir -p "${MMA_PATH}" "${MULTISCALE_WORK_DIR}" \
    && chown -R "${USERNAME}:${USERNAME}" "${MMA_PATH}" "${MULTISCALE_WORK_DIR}"

# Install container helper commands.
COPY docker/entrypoint.sh /usr/local/bin/mma-entrypoint
COPY docker/mma-build /usr/local/bin/mma-build
COPY docker/mma-jupyter /usr/local/bin/mma-jupyter
RUN chmod +x \
        /usr/local/bin/mma-entrypoint \
        /usr/local/bin/mma-build \
        /usr/local/bin/mma-jupyter

USER ${USERNAME}
WORKDIR ${MMA_PATH}

ENTRYPOINT ["mma-entrypoint"]
CMD ["bash"]
