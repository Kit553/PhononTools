#!/usr/bin/env bash

# PALMA CPU/MPI starting profile; benchmark and resize for your actual system.
# See ../../../../00_Common/00_SettingUp/VASP/benchmarking.md. Memory is per node, not per rank.
# Exclusive allocation reserves the node even when this job uses fewer ranks.
#SBATCH --job-name=gamma-dfpt
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --partition=normal
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --no-requeue
#SBATCH --account=uni
#SBATCH --mail-type=NONE
# To enable email, set --mail-type=FAIL,END and add your own --mail-user.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the calculation directory}"
ulimit -s unlimited

# Existing VASP 6 example; confirm availability and use a compatible module stack.
# The documented 2019a / VASP/5.4.4-vtst stack is a separate benchmark option.
module load palma/2020b
module load intel/2020b
module load VASP/6.2.0

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

for input in INCAR POSCAR POTCAR KPOINTS; do
    [[ -s "$input" ]] || { echo "Missing or empty input: $input" >&2; exit 1; }
done

command -v vasp_std >/dev/null
module list
srun vasp_std > vasp.out 2> vasp.err
