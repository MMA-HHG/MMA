### This shell scripts sets all the necessary env variables

### To run this script make sure you made the script executable using ```chmod```
### command and run with the command ```source ./set_env_vars.sh``` or simply as
### ```. ./set_env_vars.sh```.

### Default path with directories - should contain all the necessary git repos
export MMA_PATH="$(git rev-parse --show-toplevel)"

export HPC=Curta

export PYTHONPATH=$PYTHONPATH:$MMA_PATH/shared_python:$MMA_PATH/CUPRAD/python

export CUPRAD_HOME=$MMA_PATH/CUPRAD

export CUPRAD_BINARY=$MMA_PATH/CUPRAD/binary

export CUPRAD_SCRIPTS=$MMA_PATH/CUPRAD/scripts

export CUPRAD_INPUTS=$MMA_PATH/CUPRAD/testing

export CUPRAD_PYTHON=$MMA_PATH/CUPRAD/python

export TDSE_1D_SOURCE=$MMA_PATH/1DTDSE/

export TDSE_1D_POST_PROCESSING=$MMA_PATH/1DTDSE/post_processing

export TDSE_1D_SCRIPTS=$MMA_PATH/1DTDSE/scripts

export TDSE_1D_SLURM=$MMA_PATH/1DTDSE/slurm

export TDSE_1D_BINARY=$MMA_PATH/1DTDSE/binary

export HANKEL_HOME=$MMA_PATH/Hankel

export TDSE_1D_HOME=$MMA_PATH/1DTDSE

export MULTISCALE_HOME=$MMA_PATH

export MULTISCALE_SCRIPTS=$MMA_PATH/multiscale/scripts

export MULTISCALE_WORK_DIR=$MMA_PATH/work_dir

if [[ $(module > /dev/null 2>&1) ]]; then
    source $MMA_PATH/Modules/load_modules.sh
fi
