# Testing 
 - This directory contains files I used for testing CUPRAD OCC on Eclipse.

## Files
### my bashrc
 - This file contains content of the .bashrc file I had on Eclipse.
 - The most important function is ptc which stands for 'please test cuprad'.
    - It is also an alias which uses file inp.txt as an input, so I do not need to use standard input to write the input values.
 - This is the content of the file:
``` bash
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

function load_modules() {
	echo "Loading necessary modules"
	module load GCC/6.1.0
	module load OpenMPI/1.10.2-GCC-6.1.0-psm-dbg
	module load HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm
	module load FFTW/3.3.5-OpenMPI-1.10.2-GCC-6.1.0-psm-dbg
	module list
}
load_modules
cd CUPRAD_OCC/sources

function ptc() { # please test cuprad
	make code
	cp test.h5 results.h5                                                                                        
	./make_start_occigen.e
	mpirun -n 2 ./cuprad_occigen.e
}

alias ptc="ptc < inp.txt"
```

### inp.txt
 - This file contains strings needed for input to pre-processor.
 - This is the content of my version of the file:
``` bash
results.h5
0
0
0
```

### makefile
 - This is my makefile that I used to compile my code.

### test.h5
 - This file is necessary for the testing as it contains just inputs group which was created by universal input python script.
