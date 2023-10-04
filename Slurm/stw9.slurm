#!/bin/bash
#SBATCH --job-name=stw_s11
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                               # Time limit hrs:min:sec
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --output=stw_s11.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/data/troxellab/danniw/r/BS

Rscript --vanilla stepped_wedge_random_walk_gaussian_process_cluster.R