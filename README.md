# Development notes
The development of the coupled model

## PYTHONPATH organisation

Several Python procedures require modules from various places. The actual approach is to include all the locations in PYTHONPATH

- There is needed to inlcude the independent git repository `python_modules`
- The path should also contain `shared_python` from the root of this repository
- There is an example of the paths added into the environment. The environment variables are used in some Python scripts, so they shall be added into `.bashrc` to ensure compatibility:
``` bash
export UNIV_INPUT_PATH=/home/vabekjan/git/universal_input

export PYTHONPATH=$PYTHONPATH:/home/vabekjan/git/python_modules:/home/vabekjan/git/CUPRAD_TDSE_Hankel/shared_python

export CUPRAD_HOME=/home/vabekjan/git/CUPRAD_TDSE_Hankel/CUPRAD

export CUPRAD_BINARY=/home/vabekjan/git/CUPRAD_TDSE_Hankel/CUPRAD/binary

export CUPRAD_SCRIPTS=/home/vabekjan/git/CUPRAD_TDSE_Hankel/CUPRAD/scripts

export CUPRAD_INPUTS=/home/vabekjan/git/CUPRAD_TDSE_Hankel/CUPRAD/testing

export CUPRAD_PYTHON=/home/vabekjan/git/CUPRAD_TDSE_Hankel/CUPRAD/python

export TDSE_1D_SOURCE=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE/

export TDSE_1D_POST_PROCESSING=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE/post_processing

export TDSE_1D_SCRIPTS=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE/scripts

export TDSE_1D_SLURM=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE/slurm

export TDSE_1D_BINARY=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE/binary

export HANKEL_HOME=/home/vabekjan/git/CUPRAD_TDSE_Hankel/Hankel

export TDSE_1D_HOME=/home/vabekjan/git/CUPRAD_TDSE_Hankel/1DTDSE

export MULTISCALE_HOME=/home/vabekjan/git/CUPRAD_TDSE_Hankel

export MULTISCALE_SCRIPTS=/home/vabekjan/git/CUPRAD_TDSE_Hankel/multiscale/scripts

export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON
```
- The part of the path `/home/vabekjan/git/` should be replaced by the paths you actually use.



## Merging independent repos into the multi-scale model:
https://stackoverflow.com/a/52933095