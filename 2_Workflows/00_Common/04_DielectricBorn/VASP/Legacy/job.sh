#!/bin/bash

#SBATCH --job-name=Dielectric
###SBATCH --output=$SLURM_SUBMIT_DIR/out_job_%j
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --partition=normal
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --no-requeue
#SBATCH --exclusive
#SBATCH --account=uni

#SBATCH --mail-type=NONE
#SBATCH --mail-user=user@uni-muenster.de # your mail address

module load palma/2020b
module load intel/2020b
module load VASP/6.2.0

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

cd $SLURM_SUBMIT_DIR
srun vasp_std > out_VASP
