#!/bin/bash
#SBATCH --job-name=stw_s45
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_dev
#SBATCH --time=8:00:00                               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB
#SBATCH --output=stw_s45
#SBATCH --array=1-150                             # 100 jobs in the array

module load r/4.2.2
cd /gpfs/data/troxellab/danniw/r/BS

Rscript --vanilla decreasing_cluster_specific_time_varying_intervention.R $SLURM_ARRAY_TASK_ID