#!/bin/bash -l
#$ -pe mpi 16
#$ -q development.q
#$ -N test_nallo
#$ -cwd
#$ -t 1-7

module load miniconda
source activate nf-core

nf-test test --tag PIPELINE --shard $SGE_TASK_ID/7 --update-snapshot --clean-snapshot --profile +singularity
