#!/bin/bash
#PBS -l select=1:ncpus=1:mem=500mb
#PBS -l walltime=00:59:00
#PBS -A "Coilseq"

## adapted from https://help.igb.illinois.edu/Job_Array_Example

set -e # exit immediately upon error

module load R # LOAD R

# set PATH VARIABLES
R_LIBS=/home/machne/R

#### TODO: generate param list in loop (as for SGE qsub)
#### and write to job.conf
## numj=$(expr `cat job.conf | wc -l`)
## qsub -o $logdir -e $logdir -J 1-${numj} scripts/pbs_batch.sh

jobList="job.conf"
job=$( head -n $PBS_ARRAYID $jobList | tail -n 1 )
echo $job
$job > $logdir/${PBS_JOBID}_${PBS_ARRAYID}.dat 2> $logdir/${PBS_JOBID}_${PBS_ARRAYID}.log


