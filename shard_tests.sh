#!/bin/bash

#$ -N shard_tests
#$ -cwd
#$ -pe mpi 16
#$ -q development.q
#$ -V
#$ -t 1-10

nf-test test --profile +singularity --shard $SGE_TASK_ID/$SGE_TASK_LAST --update-snapshot --tag PIPELINE
