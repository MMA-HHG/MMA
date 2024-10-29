####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2020)
#
# The purpose of this code is to merge outputs generated from 1DTDSE. It uses a general
# procedure for assigning indices based on keys, it thus can be easily generalised for more
# dimensions.



# finding mn from external location
# https://stackoverflow.com/questions/17198319/how-to-configure-custom-pythonpath-with-vm-and-pycharm

import sys
import h5py
import MMA_administration as MMA



arguments = sys.argv
filename = arguments[1]
with h5py.File(filename,'a') as f:
    del f[MMA.paths['CTDSE_outputs']]   
print('Done')