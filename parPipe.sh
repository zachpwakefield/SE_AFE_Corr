#!/bin/bash -l


#$ -P evolution
#$ -l h_rt=99:59:00
#$ -N parPipeKend
#$ -j y
#$ -m ea
#$ -pe omp 28
#$ -l mem_per_core=18G

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load bedops bedtools

module load R/4.1.1


##change to whatever is necessary
# sample_directory = /projectnb2/evolution/zwakefield/finPipe/runs/
# output_directory = /projectnb2/evolution/zwakefield/finPipe/Corr/corrOut/

Rscript /projectnb2/evolution/zwakefield/finPipe/Corr/parPipeKend.R /projectnb2/evolution/zwakefield/finPipe/runs/ /projectnb2/evolution/zwakefield/finPipe/Corr/corrOut/
