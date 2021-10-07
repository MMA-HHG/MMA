#!/bin/bash

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

module purge
module load intel/17.2 python/3.6.3
python3 $UNIV_INPUT_PATH/process_multiparametric.py -i-reg $regular -i-mp $multiparametric -multiparam-groups -keep-intermediate -univ-inps -ohdf5 results.h5 -g inputs

cp multiparameters/*.h5 .


for ((iter=init; iter <= final ; iter++))
do

    simulation=results_${iter}.h5
    ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
    mkdir sim_$ksimulation
    mv $simulation sim_$ksimulation/
    cd sim_$ksimulation
        echo "executing simulation $ksimulation"
        $CUPRAD_SCRIPTS/run_CUPRAD.sh
    cd ..

done


# for ((iter=init; iter <= final ; iter++))
# do
#    simulation=results_${iter}.h5
#    ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
#    mkdir sim_$ksimulation
#    cd sim_$ksimulation
#         echo "executing simulation $ksimulation"
#         echo "h5name $simulation"
#    cd ..

# done

# for simulation in results_*.h5; do

#     ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
#     mkdir sim_$ksimulation
#     mv $simulation sim_$ksimulation/
#     cd sim_$ksimulation
#         echo "executing simulation $ksimulation"
#         $CUPRAD_SCRIPTS/run_CUPRAD.sh
#     cd ..

# done;


echo "done" 
