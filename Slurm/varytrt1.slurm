#!/bin/bash
#SBATCH --job-name=stw_s30
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=24:00:00                               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --output=stw_s30.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/home/dw2625/r/BS

Rscript --vanilla stepped_wedge_time_vary_trt_m1.R