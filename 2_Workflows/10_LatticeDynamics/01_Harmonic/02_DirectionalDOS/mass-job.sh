#!/usr/bin/env bash

#SBATCH --job-name={{JOB_NAME}}
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G

set -euo pipefail

# Adapt resources and activate the CAT/Phonopy/mpi4py environment here.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd {{WORKDIR}}

srun {{FREQBALL_COMMAND}}
