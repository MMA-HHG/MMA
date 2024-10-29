#!/bin/bash

# # rm -r * # Clean folder # rather do it manually

while getopts :r:m:i:f: flags; do
  case "${flags}" in
    i)
      init=${OPTARG}
      echo "i: $init"
      ;;
    f)
      final=${OPTARG}
      echo "f: $final"
      ;;
    r)
      regular=${OPTARG}
      echo "regular input: $regular"
      ;;
    m)
      multiparametric=${OPTARG}
      echo "multiparametric input: $multiparametric"
      ;;
    :)
      echo "Option -$OPTARG requires an argument."
    #   exit 1
      ;;
  esac
done
# echo "test" 

for i in {1..5}
do
   echo "Welcome $i times"
done

# init=5
# final=7

for ((a=init; a <= final ; a++))
do
   echo "Welcome $a times"
done



# for i in {$((init))..$((final))}
# do
#    echo "Welcome $i times"
# done

# exit 1

# cp $TESTPATH/*.inp .

# rm -r sim_*
# rm *.h5

# module purge
# module load intel/17.2 python/3.6.3
# python3 $UNIV_INPUT_PATH/process_multiparametric.py -i-reg $regular -i-mp $multiparametric -multiparam-groups -keep-intermediate -univ-inps -ohdf5 results.h5 -g inputs

# cp multiparameters/*.h5 .


# for simulation in results_*.h5; do

    # ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
    # mkdir sim_$ksimulation
    # mv $simulation sim_$ksimulation/
    # cd sim_$ksimulation
    #     echo "executing simulation $ksimulation"
    #     # $CUPRAD_SCRIPTS/run_CUPRAD.sh
    # cd ..

# done;

echo "done" 
