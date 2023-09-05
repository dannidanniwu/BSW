#!/bin/bash
#SBATCH --job-name=stw_s3
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_medium
#SBATCH --time=24:00:00                               # Time limit hrs:min:sec
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --output=stw_s3.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/data/troxellab/danniw/r

Rscript --vanilla GAM_timespline_bayes_vs_freq.R