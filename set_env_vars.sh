### This shell scripts sets all the necessary env variables

### To run this script make sure you made the script executable using ```chmod``` 
### command and run with the command ```source ./set_env_vars.sh``` or simply as
### ```. ./set_env_vars.sh```.

### Default path with directories - should contain all the necessary git repos
export GIT_PATH=/gpfs/home/tnemec

export HPC=Curta

export UNIV_INPUT_PATH=$GIT_PATH/universal_input

export MSM_PATH=$GIT_PATH/CUPRAD_TDSE_Hankel

export PYTHONPATH=$PYTHONPATH:$GIT_PATH/python_modules:$MSM_PATH/shared_python

export CUPRAD_HOME=$MSM_PATH/CUPRAD

export CUPRAD_BINARY=$MSM_PATH/CUPRAD/binary

export CUPRAD_SCRIPTS=$MSM_PATH/CUPRAD/scripts

export CUPRAD_INPUTS=$MSM_PATH/CUPRAD/testing

export CUPRAD_PYTHON=$MSM_PATH/CUPRAD/python

export TDSE_1D_SOURCE=$MSM_PATH/1DTDSE/

export TDSE_1D_POST_PROCESSING=$MSM_PATH/1DTDSE/post_processing

export TDSE_1D_SCRIPTS=$MSM_PATH/1DTDSE/scripts

export TDSE_1D_SLURM=$MSM_PATH/1DTDSE/slurm

export TDSE_1D_BINARY=$MSM_PATH/1DTDSE/binary

export HANKEL_HOME=$MSM_PATH/Hankel

export TDSE_1D_HOME=$MSM_PATH/1DTDSE

export MULTISCALE_HOME=$MSM_PATH

export MULTISCALE_SCRIPTS=$MSM_PATH/multiscale/scripts

source Modules/load_modules.sh
