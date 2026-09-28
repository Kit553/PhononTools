#!/bin/bash

#SBATCH --job-name=Na3PS4_tet_PhonoPy         # the name of your job
###SBATCH --output=$SLURM_SUBMIT_DIR/out_job_%j
#SBATCH --nodes=1                # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=32        # the number of CPU cores per node
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=08:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --mem-per-cpu=2G         #RAM per CPU
#SBATCH --no-requeue
#SBATCH --exclusive
#SBATCH --account=uni

#SBATCH --mail-type=NONE             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=peckert@uni-muenster.de # your mail address

module load palma/2020b
module load intel/2020b
module load VASP/6.2.0

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

cd $SLURM_SUBMIT_DIR
srun vasp_std > out_VASP
