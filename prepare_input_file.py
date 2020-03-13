####################################################################################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE) (2020)
#
# the purpose of this code is to generate a universal input HDF5-file,
# it allows to write all the inputs in a free form specified by
#
#	'name'	'value' 'type'	'unit'	'comments, etc. (no effect on the file)'
#
# type is in (R)eal (I)nteger (S)tring
#
# The idea is to create universal input files that can be used accross codes. 
# All the data are now stored in the group 'inputs'. Be carefull with the names of
# variables. If they are chosen correctly, it would be easier to cross-platform it.

import numpy as np
import h5py

#### THE MAIN PROGRAM #####


inputfilename = 'FreeFormInputs.inp'
outputfilename = 'results.h5'


## specify the name of files and datasets
InputFile = open(inputfilename,"r")
GeneratedFile = h5py.File(outputfilename, 'w') # be careful, now it trucates if already exists
grp = GeneratedFile.create_group('inputs')

def adddataset(h_path,sep_line,line):
  if ( (0 < len(sep_line)) and (len(sep_line) < 3) ): print('warning in the line (input too short, entry skipped): ' + line); return

  if (sep_line[2]=='R'): dset_id = h_path.create_dataset(sep_line[0], data=float(sep_line[1]))
  elif (sep_line[2]=='I'): dset_id = h_path.create_dataset(sep_line[0], data=int(sep_line[1]))
  elif (sep_line[2]=='S'): dset_id = h_path.create_dataset(sep_line[0], data=np.string_(sep_line[1]))
  else: print('warning in the line (type unrecognised and entry skipped): ' + line); return

  if (len(sep_line) < 4):
    print('warning in the line (missing units ?): ' + line) 
    dset_id.attrs['units']=np.string_('?') 
  else: dset_id.attrs['units']=np.string_('['+sep_line[3]+']')



lines = InputFile.readlines()
for line in lines:
  sep_line = line.split(); # separate the line
  if ( (len(sep_line)==0) or (sep_line[0]=='#') or (sep_line[0]=='##') ): pass #print('empty or commented line')
  else: adddataset(grp,sep_line,line)


InputFile.close()
GeneratedFile.close()


print('done');
