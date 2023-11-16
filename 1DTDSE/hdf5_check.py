"""
This python script checks if the HDF5 file contains necessary inputs for the 
1DTDSE.
"""

import h5py
import argparse

def test_key(key):
    try:
        f[key][()]
    except KeyError:
        print("Missing key: ", key)
        return key
    else:
        pass


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--inhdf5", required=False, help="Path to HDF5 file to check.")
args = vars(ap.parse_args())

filename = None

if args['inhdf5'] == None:
    print("Type the HDF5 file to check with relative path: ")
    filename = input()
else:
    filename = args['inhdf5']

with h5py.File(filename, 'r') as f:
    keys = ["TDSE_inputs/Eguess",
              "TDSE_inputs/N_r_grid",
              "TDSE_inputs/dx",
              "TDSE_inputs/InterpByDTorNT",
              "TDSE_inputs/dt",
              "TDSE_inputs/Ntinterp",
              "TDSE_inputs/analy_writewft",
              "TDSE_inputs/analy_tprint",
              "TDSE_inputs/x_int",
              "TDSE_inputs/trg_a",
              "TDSE_inputs/CV_criterion_of_GS",
              "TDSE_inputs/gauge_type",
              ]

    missing_keys = []

    for key in keys:
        out = test_key(key)
        if out is not None:
            missing_keys.append(out)

    print("*******************************\n")
    if missing_keys!=[]:
        print("Missing keys:")
        print(missing_keys)
    else:
        print("All OK.")
    print("Check finished.")