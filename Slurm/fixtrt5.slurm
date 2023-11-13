#!/bin/bash
#SBATCH --job-name=stw_s28
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=24:00:00                               # Time limit hrs:min:sec
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=28GB
#SBATCH --output=stw_s28.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/home/dw2625/r/BS

Rscript --vanilla stepped_wedge_random_walk_fixed_trt1_local.R