#!/bin/bash



test_presets=("test1_modT1" "test1_mod_incT1" "test1_mod_decT1"
              "test1_modT2" "test1_mod_incT2" "test1_mod_decT2"
              "test1_modT3" "test1_mod_incT3" "test1_mod_decT3"
              "test2_modT1" "test2_mod_incT1" "test2_mod_decT1"
              "test2_modT2" "test2_mod_incT2" "test2_mod_decT2"
              "test2_modT3" "test2_mod_incT3" "test2_mod_decT3")



echo ${test_presets[1]}
echo ${test_presets[4]}



module purge
module load hdf5
module load intel-mkl

for k1 in ${!test_presets[@]}; do

echo "----------------------------------------------------------------------"
echo ${test_presets[k1]}
echo "----------------------------------------------------------------------"

mkdir ${test_presets[k1]}
cd ${test_presets[k1]}
$CUPRAD_HOME/build/make_start.e <<INPUTS
${test_presets[k1]}
0
0
0
INPUTS

sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm

cd ..

done



mkdir test1_GaussT2

cd test1_GaussT2
$CUPRAD_HOME/build/make_start.e <<INPUTS
${test_presets[3]}
0
0
0
INPUTS

module purge
module load python
python3 $CUPRAD_HOME/benchmarks/add_dens_modulation1.py 1e-5

cd ..


mkdir test2_GaussT2

cd test2_GaussT2

module purge
module load hdf5
module load intel-mkl

$CUPRAD_HOME/build/make_start.e <<INPUTS
${test_presets[12]}
0
0
0
INPUTS

module purge
module load python
python3 $CUPRAD_HOME/benchmarks/add_dens_modulation1.py 5e-4
sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm
cd ..


echo "done" 
