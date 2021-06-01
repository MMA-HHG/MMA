import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
import mynumerics as mn

# import mynumerics as mn
import matplotlib.pyplot as plt


#
arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    results_path = os.path.join("D:\data", "Discharges")
    results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "foci")
    results_path = os.path.join("D:\data", "Discharges", "f_scan")
    results_path = os.path.join("D:\data", "Discharges", "TDSE", "t6")
    
    results_TDSE = os.path.join("D:\data", "Discharges", "TDSE", "TDSE10planes1")


file = 'results_1.h5' # 'results_Ar_vac.h5', 'Ar_vac_long.h5' 'results_3.h5' 'results_1.h5'

file_TDSE = 'results_merged.h5' # 'hdf5_temp_0000000.h5'

file_TDSE = os.path.join(results_TDSE,file_TDSE)

index_select = [(0,0),(10,0),(0,2)] # (r,z)
kz_TDSE = 0
kr_TDSE = 0
kz_CUPRAD = 0
kr_CUPRAD = 0

r_select = [k1[0] for k1 in index_select]
z_select = [k1[1] for k1 in index_select]


def get_slices(dset,i_select):
    N_select = len(i_select)
    dset_shape = dset.shape
    if (len(dset_shape) == 3):
        res = np.zeros((N_select,dset_shape[-1]))
        for k1 in range(N_select):
            res[k1,:] = dset [i_select[k1][0],i_select[k1][1],:]
    if (len(dset_shape) == 4):  
        res = np.zeros((N_select,dset_shape[-1]), dtype = np.cdouble)
        for k1 in range(N_select):
            res[k1,:] = dset [i_select[k1][0],i_select[k1][1],:,0] + \
                        1j*dset [i_select[k1][0],i_select[k1][1],:,1]
    return res

file_path = os.path.join(results_path,file)
print('processing:', file_path)             
with h5py.File(file_TDSE, 'r') as InputArchiveTDSE:
   
   Efield_full =  InputArchiveTDSE['Efield'][:]
   
   omega0 = 0.05752948549410085
   # Efield_TDSE = InputArchiveTDSE['Efield'][r_select,z_select,:]
   
   Efield_TDSE = Efield_full[r_select,z_select,:]
   
   field_shape = InputArchiveTDSE['Efield'].shape
   
   Efield_selected = get_slices(InputArchiveTDSE['Efield'],index_select)
   
   # FEfield_TDSE = InputArchiveTDSE['FEfield'][kr_TDSE,kz_TDSE,:,0] + \
   #                     1j*InputArchiveTDSE['FEfield'][kr_TDSE,kz_TDSE,:,1]
   # tgrid_TDSE = InputArchiveTDSE['tgrid'][:]
   
   # SourceTerm_TDSE = InputArchiveTDSE['SourceTerm'][kr_TDSE,kz_TDSE,:]
   # FSourceTerm_TDSE = InputArchiveTDSE['FSourceTerm'][kr_TDSE,kz_TDSE,:,0] + \
   #                    1j*InputArchiveTDSE['FSourceTerm'][kr_TDSE,kz_TDSE,:,1]
   # ogrid_TDSE = InputArchiveTDSE['omegagrid'][:]
   # PopTot_TDSE = InputArchiveTDSE['PopTot'][kr_TDSE,kz_TDSE,:]
   # PopInt_TDSE = InputArchiveTDSE['PopInt'][kr_TDSE,kz_TDSE,:]
   # expval_x_TDSE = InputArchiveTDSE['expval_x'][kr_TDSE,kz_TDSE,:]
   
   # GS_init = InputArchiveTDSE['ground_state'][:,0] + 1j*InputArchiveTDSE['ground_state'][:,1]
   # xgrid_micro = InputArchiveTDSE['xgrid_micro'][:]




# fig = plt.figure()
# plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE)
# plt.title('TDSE')
# plt.show()
# # plt.close(fig)

# # fig = plt.figure()
# # plt.plot(tgrid_CUPRAD-tgrid_CUPRAD[0],Efield_CUPRAD)
# # # plt.plot(tgrid_TDSE*units.TIMEau,Efield_TDSE*units.EFIELDau)
# # plt.title('TDSE, CUPRAD')
# # plt.show()
# # # plt.close(fig)

# fig = plt.figure()
# # plt.plot(tgrid_TDSE,Efield_TDSE)
# plt.plot(tgrid_TDSE,SourceTerm_TDSE)
# plt.title('Field, SourceTerm')
# plt.show()
# # plt.close(fig)

# fig = plt.figure()
# plt.plot(tgrid_TDSE,
#          mn.apply_filter(SourceTerm_TDSE, mn.filter_box, tgrid_TDSE, [6750,tgrid_TDSE[-1]])
#          )
# plt.title('SourceTerm, filtered')
# plt.show()
# # plt.close(fig)

# fig = plt.figure()
# plt.plot(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
# plt.title('FSourceTerm')
# plt.show()
# # plt.close(fig)

# fig = plt.figure()
# plt.semilogy(ogrid_TDSE/omega0,abs(FSourceTerm_TDSE))
# plt.title('FSourceTerm')
# plt.show()
# # plt.close(fig)

# # filter spectrum
# ogrid_FE, FE_filter, Nt = mn.fft_t_nonorm(
#                             tgrid_TDSE,
#                             mn.apply_filter(SourceTerm_TDSE, mn.filter_box, tgrid_TDSE, [6750,tgrid_TDSE[-1]])
#                             )

# fig = plt.figure()
# plt.semilogy(ogrid_FE/omega0,abs(FE_filter))
# plt.title('FS_filt')
# plt.show()
# # plt.close(fig)

# ogrid_FE, FE_filter, Nt = mn.fft_t_nonorm(
#                             tgrid_TDSE,
#                             mn.apply_filter(SourceTerm_TDSE, mn.filter_box, tgrid_TDSE, [6750,tgrid_TDSE[-1]], apply_blackman = True)
#                             )

# fig = plt.figure()
# plt.semilogy(ogrid_FE/omega0,abs(FE_filter))
# plt.title('FS_filt, blackman')
# plt.show()
# # plt.close(fig)


# fig = plt.figure()
# # plt.semilogy(ogrid_TDSE/omega0,abs(FEfield_TDSE))
# plt.title('FEfield')
# plt.show()
# # plt.close(fig)

# fig = plt.figure()
# plt.plot(tgrid_TDSE,PopTot_TDSE)
# plt.title('PopTot')
# plt.show()

# fig = plt.figure()
# plt.plot(tgrid_TDSE,PopInt_TDSE)
# plt.title('PopInt')
# plt.show()

# fig = plt.figure()
# plt.plot(tgrid_TDSE,expval_x_TDSE)
# plt.title('<x>')
# plt.show()
# # plt.close(fig)

# # fig = plt.figure()
# # plt.semilogy(tgrid_TDSE,PopTot_TDSE)
# # plt.title('PopTot')
# # plt.show()
# # # plt.close(fig)

# fig = plt.figure()
# plt.plot(xgrid_micro,abs(GS_init))
# plt.title('GS')
# plt.show()
# # plt.close(fig)

# fig = plt.figure()
# plt.semilogy(xgrid_micro,abs(GS_init))
# plt.title('GS')
# plt.show()
# # plt.close(fig)

# # Gabor
# # t_G, o_G, Gabor = mn.gabor_transf(SourceTerm_TDSE, tgrid_TDSE, 5000, 10000, 2000, 8.0)
# t_G, o_G, Gabor = mn.gabor_transf(SourceTerm_TDSE, tgrid_TDSE, 7000, 9000, 2000, 8.0, omegamax=40*omega0)

# o_Gr, t_Gr, Gaborr = mn.interpolate_2D(o_G,t_G,Gabor,2000,2000)

# vmin = np.max(np.log(Gaborr))-6.
# fig = plt.figure()
# plt.pcolor(t_Gr,o_Gr/omega0,(np.log(Gaborr)).T, shading='auto',vmin=vmin)
# plt.title('Gabor')
# plt.show()
# # plt.close(fig)



