#!/bin/bash
#SBATCH --job-name=LUMA_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --time=0:10:0
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=ecseaa28

source ./LUMA/scripts/ARCHER2/setup.sh

export OMP_NUM_THREADS=16

# srun to launch the executable
srun --cpu-bind=cores ./LUMA > luma.log
