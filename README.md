# Maker fringes  project

## Parameters

- argon
    - $I = 1 \times 10^{14}~\mathrm{W/cm^2}$
    - $w_0 = 25~\mathrm{\mu m}$
    - $\lambda = 800~\mathrm{nm}$
    - $\tau_{\text{FWHM}} = 30~\mathrm{fs}$
    - Gaussian mode (TEM00)
    - $p = 0.1~\mathrm{bar}$
    - $L_{\text{med}} = 50~\mathrm{\mu m}$
    - focusing: $\pm 100~\mathrm{\mu m}$ w.r.t. the focus
    - $z_{FF} = 50~\mathrm{cm}$

- xenon
    - $I = 7 \times 10^{13}~\mathrm{W/cm^2}$
    - $w_0 = 10~\mathrm{\mu m}$
    - $\lambda = 1030~\mathrm{nm}$
    - $\tau_{\text{FWHM}} = 200~\mathrm{fs}$
    - Gaussian mode (TEM00)
    - $p = 0.1~\mathrm{bar}$
    - $L_{\text{med}} = 50~\mathrm{\mu m}$
    - focusing: $\pm 100~\mathrm{\mu m}$ w.r.t. the focus
    - $z_{FF} = 50~\mathrm{cm}$

<!-- - xenon
30 fs - 800 nm
200 fs - 1030 nm

Conditions are:
No envelop temporally
Integration on 1 optical cycle
Trajectory selection (Lewenstein model)
I0=1e14 W/cm2 // future 7e13
w0=25 microns // future 10
Lambda0=800 nm // future 1030 nm
Spatial mode= TEM00 gaussian
Pressure=0.1 bar
Atom: argon// future xenon
Absorption=1cm response
Length of medium= 50 micron
Position of target= @ focus and before z=-100 micron and after z=+100
Calculation plan: @ farfield this 50 cm from medium output -->


# Development notes
The development of the coupled model

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