#!/bin/bash



test17_presets=( "test17_modT1" "test17_modT1" "test17_modT1")
test18_presets=("test18_modT1" "test18_modT1" "test18_modT1")
test_1_jet=("jet1_GaussT1" "jet1_GaussT2" "jet1_GaussT3")
test_2_jet=("jet2_GaussT1" "jet2_GaussT2" "jet2_GaussT3")
test_10_jet=("jet10_GaussT1" "jet10_GaussT2" "jet10_GaussT3")


# 1 jet
for k1 in ${!test_1_jet[@]}; do

echo "----------------------------------------------------------------------"
echo ${test_1_jet[k1]}
echo "----------------------------------------------------------------------"

mkdir ${test_1_jet[k1]}
cd ${test_1_jet[k1]}

module purge
module load hdf5
module load intel-mkl

$CUPRAD_HOME/build/make_start.e <<INPUTS
${test17_presets[k1]}
0
0
0
INPUTS

module purge
module load python/3.9.12
python3 $CUPRAD_HOME/benchmarks/add_dens_modulation2.py

sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm

cd ..

done


# 2 jet
for k1 in ${!test_2_jet[@]}; do

echo "----------------------------------------------------------------------"
echo ${test_2_jet[k1]}
echo "----------------------------------------------------------------------"

mkdir ${test_2_jet[k1]}
cd ${test_2_jet[k1]}

module purge
module load hdf5
module load intel-mkl

$CUPRAD_HOME/build/make_start.e <<INPUTS
${test17_presets[k1]}
0
0
0
INPUTS

module purge
module load python/3.9.12
python3 $CUPRAD_HOME/benchmarks/add_dens_modulation3.py

sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm

cd ..

done



# 10 jet
for k1 in ${!test_10_jet[@]}; do

echo "----------------------------------------------------------------------"
echo ${test_10_jet[k1]}
echo "----------------------------------------------------------------------"

mkdir ${test_10_jet[k1]}
cd ${test_10_jet[k1]}

module purge
module load hdf5
module load intel-mkl

$CUPRAD_HOME/build/make_start.e <<INPUTS
${test18_presets[k1]}
0
0
0
INPUTS

module purge
module load python/3.9.12
python3 $CUPRAD_HOME/benchmarks/add_dens_modulation10.py

sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm

cd ..

done


echo "done" 
