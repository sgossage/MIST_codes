#!/bin/bash

#SBATCH -n <<CORES>>
#SBATCH -N 1
#SBATCH -t <<RUNTIME>>
#SBATCH --mem 8000
#SBATCH -p <<PARTITION>>
#SBATCH -o <<RUNNAME>>.o
#SBATCH -e <<RUNNAME>>.e

export OMP_NUM_THREADS=8
cd <<DIRNAME>>
./clean
./mk
echo "SLURM JOB ID: $SLURM_JOB_ID"
./rn
