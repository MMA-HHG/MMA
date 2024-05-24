import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
import units
import mynumerics as mn

# import mynumerics as mn
import matplotlib.pyplot as plt




arguments = sys.argv

try:
    filename = arguments[1]
    kr = int(arguments[2])
    kz = int(arguments[3])
except:
    filename = 'results_merged_example.h5'
    kr = 0
    kz = 0 



showplots = not('-nodisplay' in arguments)


# sys.exit()

print('processing:', filename)             
with h5py.File(filename, 'r') as InputArchive:
    # load data
   
   Efield = InputArchive['/Efield'][kr,kz,:];


print('data loaded:')

fig10, ax10 = plt.subplots()
ax10.plot(Efield)

# ax10.set_xlabel('r [mum]'); ax10.set_ylabel('I [cutoff]'); ax10.set_title('t=0 fs, intensity'+title_string)   
fig10.savefig('Efield_plot.png', dpi = 600)           

if showplots: plt.show()
else: plt.close()






