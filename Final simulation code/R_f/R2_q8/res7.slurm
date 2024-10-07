#!/bin/bash
#SBATCH --job-name=res_7
#SBATCH --mail-type=END,FAIL                         # send email if the job end or fail
#SBATCH --mail-user=danni.wu@nyulangone.org
#SBATCH --partition=cpu_dev
#SBATCH --time=12:00:00                               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=2
#SBATCH --mem=4GB
#SBATCH --output=res_7
#SBATCH --array=1-150                             # 100 jobs in the array

module load r/4.2.2
cd /gpfs/home/dw2625/r/Review_sim/rcode

Rscript --vanilla inc_dec_time_vary_trt_m2_pointwise_ci.R $SLURM_ARRAY_TASK_ID