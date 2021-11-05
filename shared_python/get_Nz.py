import numpy as np
import h5py


filename = 'results.h5'

with h5py.File(filename, 'r') as InputArchive:
    print('Nz = ', len(InputArchive['/outputs/zgrid'][:]))

print('Done')