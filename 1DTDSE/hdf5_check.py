import h5py

def test_key(key):
    try:
        f[key][()]
    except KeyError:
        print("Missing key: ", key)
        return key
    else:
        pass

with h5py.File("results.h5", 'r') as f:
    keys = ["TDSE_inputs/Eguess",
              "TDSE_inputs/N_r_grid",
              "TDSE_inputs/N_r_grid_exp",
              "TDSE_inputs/dx",
              "TDSE_inputs/InterpByDTorNT",
              "TDSE_inputs/dt",
              "TDSE_inputs/Ntinterp",
              "TDSE_inputs/textend",
              "TDSE_inputs/analy_writewft",
              "TDSE_inputs/analy_tprint",
              "TDSE_inputs/x_int",
              "TDSE_inputs/textend",
              "TDSE_inputs/PrintGaborAndSpectrum",
              "TDSE_inputs/a_Gabor",
              "TDSE_inputs/omegaMaxGabor",
              "TDSE_inputs/dtGabor",
              "TDSE_inputs/tmin1window",
              "TDSE_inputs/tmax1window",
              "TDSE_inputs/tmin2window",
              "TDSE_inputs/tmax2window",
              "TDSE_inputs/PrintOutputMethod",
              "TDSE_inputs/trg_a",
              "TDSE_inputs/CV_criterion_of_GS",
              "TDSE_inputs/gauge_type",
              ]

    missing_keys = []

    for key in keys:
        out = test_key(key)
        if out is not None:
            missing_keys.append(out)

    print("Missing keys:")
    print(missing_keys)
    print("Check finished.")