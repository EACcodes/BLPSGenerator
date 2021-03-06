#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -t 01:00:00
#SBATCH --requeue

echo "===================================================="
echo " submitted to local node: " $HOSTNAME
echo "===================================================="

ABINIT=/directory/abinit
ABINIT_BLPS=/directory/abinit-blps

cd $SLURM_SUBMIT_DIR

srun -n 8 $ABINIT_BLPS < INPUT > OUTPUT 

