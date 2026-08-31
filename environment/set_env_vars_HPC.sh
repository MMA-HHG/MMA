# This shell script sets the environment dor the MMA package.
# The environment variables needed in user's .bashrc/.bashrc/...
# are:
# MMA_PATH - the directory with the model
# HPC - the name of the computaitonal cluster with pre-defined modules etc.
# MULTISCALE_WORK_DIR - the directory for small experiments with the code:
#        (preparing inputs etc. shall be located on scratch or home,
#         usually home contains only small data! check HPC's policy)
#
# This environment is set by `source $MMA_PATH/environment/set_env_vars_HPC.sh`


# export MMA_PATH="$(git rev-parse --show-toplevel)"

# export HPC=Curta

export PYTHONPATH=$PYTHONPATH:$MMA_PATH/shared_python

export CUPRAD_HOME=$MMA_PATH/CUPRAD
export CUPRAD_BUILD=$MMA_PATH/CUPRAD/build
export CUPRAD_SCRIPTS=$MMA_PATH/CUPRAD/scripts
export CUPRAD_PYTHON=$MMA_PATH/CUPRAD/python
export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON

export TDSE_1D_HOME=$MMA_PATH/1DTDSE
export TDSE_1D_PYTHON=$MMA_PATH/1DTDSE/python
export TDSE_1D_SCRIPTS=$MMA_PATH/1DTDSE/scripts
export TDSE_1D_SLURM=$MMA_PATH/1DTDSE/slurm
export TDSE_1D_BUILD=$MMA_PATH/1DTDSE/build
export PYTHONPATH=$PYTHONPATH:$TDSE_1D_HOME

export HANKEL_HOME=$MMA_PATH/Hankel

export MULTISCALE_HOME=$MMA_PATH
export MULTISCALE_SCRIPTS=$MMA_PATH/multiscale/scripts

export FSPA_PATH=$MMA_PATH/FSPA

# export MULTISCALE_WORK_DIR=/mnt/d/data/work_dir

source $MMA_PATH/Modules/load_modules.sh
