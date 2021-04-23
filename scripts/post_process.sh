#!/bin/bash

while getopts :t: flags; do
  case $flags in
    t)
      analysis=${OPTARG}
      echo "-a was triggered, Parameter: $OPTARG"
      ;;
    # \?)
    #   echo "Invalid option: -$OPTARG"
    #   exit 1
    #   ;;
    :)
      echo "Option -$OPTARG requires an argument."
    #   exit 1
      ;;
  esac
done

# echo $analysis

case $analysis in
    general)
      sbatch $TESTPATH/coherence_map.slurm
      echo 'general submitted'
    #   echo 'plot 1'
      ;;
    compare)
      sbatch $TESTPATH/compare_n_results_XUVshift.slurm
      echo 'comparison submitted'
      ;;
esac

# sbatch $TESTPATH/plot1.slurm

echo "done" 
