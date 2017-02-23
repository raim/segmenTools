#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=01:59:00
#PBS -A "Coilseq"
 
set -e # exit immediately upon error

module load R

R_LIBS=/home/machne/R
YEASTSEQ=/home/machne/data/yeast/RNAseq/results
YEASTDAT=/home/machne/data/yeast
GENDAT=/home/machne/data/genomeData/
GENBRO=/home/machne/programs/genomeBrowser/

chrfile=$YEASTDAT/chromosomes/sequenceIndex_R64-1-1_20110208.csv
OUTSEQ=$YEASTSEQ
rundate=20170223 # primseg v5 - calculated with j<=i 
rundir=$OUTSEQ/segmentation/$rundate

segs=4441

~/programs/segmenTools/scripts/runSegmentier.R  --segs $segs -i $OUTSEQ/segmentation/primarysegments_v5/allprimseg.csv --primdir $OUTSEQ/segmentation/primarysegments_v5 --primfiles primseg_  --chrfile $chrfile --trafo ash --dc.trafo ash --dft.range 1,2,3,4,5,6,7 --K 16 --nui.thresh 0.6 --nui.cr 1:5 --scores icor --scales 1:5 --M 100 --Mn 100 --fuse.thresh 0.2 --short.name --idsuffix test --plot --fig.type pdf -o $rundir --genbro $GENBRO --gendat $GENDAT
