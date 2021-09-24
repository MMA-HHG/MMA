mkdir temp 
mv hdf5_temp_000* temp/

module purge
module load intelpython3/2019.3

python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $HANKEL_HOME/process_files/Hankel_Inputs_mp_large_all_cummulative_noabs.inp -ohdf5 inputs_Hankel_all_cummulative_noabs.h5 -g inputs
python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $HANKEL_HOME/process_files/Hankel_Inputs_mp_large_all_cummulative.inp -ohdf5 inputs_Hankel_all_cummulative.h5 -g inputs

sbatch $HANKEL_HOME/process_files/Hankel_mp_all_cummulative_noabs.slurm
sbatch $HANKEL_HOME/process_files/Hankel_mp_all_cummulative.slurm