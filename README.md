# To do list
* Code + documentation
    * *Make compilation smooth and easy on various machines (ongoing)*
    * ***CUPRAD*** (Jan + Stefan), **[Doxygen]((https://www.doxygen.nl/)) documentation**
        * Jan - `output.f90` is missing
        * Stefan - the core routines
        * **pre-processor**: potentially very tedious, do up to some extent
    * ***1D-TDSE*** (Tadeáš) **Doxygen documentation** ✔️
    * ***interactive TDSE*** (Tadeáš)
        * Ensure compatibility of the dll library
    * ***Hankel*** (Jan) **Python docstrings** (process by [sphynx](https://www.sphinx-doc.org/en/master/tutorial/index.html))
        * *refactoring*
            * Implement density modulation
            * hard-code preset gases with their scattering factors
            * data processing (either load all $\mathscr{F}$-sources or chunk reading from the drive)
    * ***Assembling the code***
        * Work with all the modules
        * Make the execution smooth (slurm dependecies etc.)
        * Fine-adjustments of all codes
            * HDF5-organisation
            * array orderings
    * ***Pre/post-processing, data plotting***
        * organise Pythonic libraries, integrare `mynumerics` int the model, ...
    * ***CodeOcean capsule*** (Tadeáš + Jan)
    * Enhance code + executables (compilation with a 'pedantic' compiler, valgrind) (to some extent)
* Prepare test cases (integrate within the CodeOcean capsule)
    * *Thin medium*
    * *Long medium* (either this or the previus with a Gaussian density-modulation.)
    * *Pre-ionised medium*
    * *Interactive python* ✔️ via Jupyter
    * *?* Some filamentation form Stefan's papers, ...
* Writing the paper (all)
    * CPC: no new physics, only the code and its usage + already published results
    * see below independent repo + Overlaf
* Test cases & tutorials
    * ***Create a YouTube video-tutorial***
        * Go through the whole code installation and execution, comment it.
    * (?) Jupyter notebook    
* Check conventions
    * The conventions in the code shall correspond the written ones in supplmentary materials. I did checks during the development, that the code works intrinsically. However, there might be a mistake in the sign of the $e^-$ charge in TDSE etc.
    * Photoelectron spectrum norm.

## Notes
Here are notes to the development of the whole model. The development of the particular modules - ***CUPRAD***, ***1D-TDSE*** and ***Hankel*** - are in the `README.md`'s in their respective directories.

The work on the publication is being done in an independet [repo](https://github.com/sskupin/XUVIR) morrored with its [Overleaf](https://www.overleaf.com/project/64f9ee5fb00b8a641fe54780).

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
