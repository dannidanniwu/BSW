#!/bin/bash
#SBATCH --job-name=stw_s5
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_short
#SBATCH --time=12:00:00                               # Time limit hrs:min:sec
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --output=stw_s5.out                        # Standard output and error log

module load r/4.1.1
cd /gpfs/home/dw2625/r/BS

Rscript --vanilla one_x_nonlinear_cluster.R