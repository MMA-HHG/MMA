# Development notes
Here are notes to the development of the whole model. The development of the particular modules - ***CUPRAD***, ***1D-TDSE*** and ***Hankel*** - are in the `README.md`'s in their respective directories.

The work on the publication is being done in an independet repo `XXX`.


## To do list
* *Make compilation smooth and easy on various machines (ongoing)*
* Clean-up the code
* Check conventions
    * The conventions in the code shall correspond the written ones in supplmentary materials. I did checks during the development, that the code works intrinsically. However, there might be a mistake in the sign of the $e^-$ charge in TDSE etc.
* Implement density modulation (CUPRAD + Hankel), implement $(1-\eta) \alpha_{\mathrm{n}} + \eta \alpha_{\mathrm{X^+}}$ contribution to dispersion from ions
* Test a full $\rho z$-pre-ionisation profile
* Organise Pythonic procedures (Hankel + various post-processings)
    * There is used the module `mynumerics`, this shall be cleaned & integrated in the multiscale model
* Enhance code + executables (compilation with a 'pedantic' compiler, valgrind)
* Test cases & tutorials
    * ***Create a YouTube video-tutorial***
        * Go through the whole code execution and comment it.
    * (?) Jupyter notebook    
* CodeOcean capsule
* (?) Documentation: the basic overview will be given in the paper, Shall we use [doxygen](https://www.doxygen.nl/) to create more technical documentation capturing the structure of the code?


# General notes
## PYTHONPATH organisation

Several Python procedures require modules from various places. The actual approach is to include all the locations in PYTHONPATH

- There is needed to inlcude the independent git repository `python_modules`
- The path should also contain `shared_python` from the root of this repository
- There is an example of the paths added into the environment. The environment variables are used in some Python scripts, so they shall be added into `.bashrc` to ensure compatibility:
``` bash
export GIT_PATH=/users/k2255939/git

# export GIT_PATH=/gpfs/home/jvabek/git

export UNIV_INPUT_PATH=$GIT_PATH/universal_input

export MSM_PATH=$GIT_PATH/CUPRAD_TDSE_Hankel

# # export TESTPATH=/home/vabekjan/git/CUPRAD_DEVELOP/testing/job_scheduler/inps

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

# # export PYTHON_CUPRAD=/home/vabekjan/git/CUPRAD_DEVELOP/python

export PYTHONPATH=$PYTHONPATH:$CUPRAD_PYTHON
```

These aliases may be separated from `.bashrc` by including them in `.bash_aliases` via
``` bash
if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi
```




## Merging independent repos into the multi-scale model:
https://stackoverflow.com/a/52933095
