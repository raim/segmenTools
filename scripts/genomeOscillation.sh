#!/bin/bash
#$ -N genomeOsc
#$ -cwd

## #$ -pe ompi 1

# run genomeOscillation.R on SGE 

## no core dump
ulimit -c 0

$TATADIR/r/genomeOscillation.R $@




# End of file
