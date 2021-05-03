#!/bin/bash

# mkdir direct
# cp file2 direct
module load Matlab
matlab -nodisplay -nodesktop -r "run post_process2.m" > figures/matlab_output

