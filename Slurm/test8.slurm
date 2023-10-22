#!/bin/bash
#SBATCH --job-name=stw_s12
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_medium
#SBATCH --time=48:00:00                               # Time limit hrs:min:sec
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --output=stw_s12.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/data/troxellab/danniw/r/BS

Rscript --vanilla stepped_wedge_random_walk_cluster.R