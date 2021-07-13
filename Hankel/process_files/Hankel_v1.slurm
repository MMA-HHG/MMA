# #!/bin/bash
# #SBATCH -J Hankel_parallel
# #SBATCH --ntasks=1 
# #SBATCH --constraint=BDW28
# #SBATCH --cpus-per-task=28
# #SBATCH --hint=multithread
# #SBATCH --time=2:00:00
# #SBATCH --output Hankel.output
# module purge
# module load intelpython3/2019.3
# srun python Hankel_Ilist_and_phenom.py

# echo "Job finished" 

#!/bin/bash
#SBATCH -J Hankel-long
#SBATCH --nodes=1
#SBATCH --constraint=[HSW24|BDW28]
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --output Hankel-%J.output

module purge
module load intelpython3/2019.3

srun python $HANKEL_HOME/Hankel_long_medium.py -here -nodisplay

echo "Job finished" 