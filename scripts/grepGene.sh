#!/usr/bin/env bash

## TODO: parse more arguments, and allow for multiple genes!

## gene: standard gene name, IME4, IME1
gene=$1 
plot=$2
range=$3
primsegs=$YEASTSEQ/segmentation/primarysegments/primseg.csv
features=/data/yeast/feature_R64-1-1_20110208_withclusters.csv
sgtype=$YEASTSEQ/segmentation/segmentTest/20161228/annotation/T.ash_D.dcash_K.16_S.icor_E.3_M.175_

seg=${sgtype}genes.csv
sas=${sgtype}antisense.csv
id=`sed 's/"//g' $features|grep gene | grep -P "\t${gene}\t"| cut -f 2`
sid=`grep $id $seg | cut -f 1|sed 's/.*\://'`
primseg=`echo $sid|sed s/_.*//`
asid=`grep -P "^$sid\t" $sas| cut -f 2`
asprimseg=`echo $asid|sed s/_.*//`
# GREP COORDINATES
primsg=`echo $primseg | sed 's/^0*//g'`
coor=`grep -P "^$primsg\t" $primsegs | cut -f 2,3,4 | sed 's/\s/,/;s/\s/:/'`
 

echo GENE $gene 
echo PRIMSEG $primseg SEGMENT $sid 
echo ANTISENSE $asprimseg SEGMENT $asid 
echo $coor

## use genomeBrowser to plot the primary segment
if  [ "$plot" = "plot" ]; then
    settings=$GENBRO/data/selections.R
    genome=$GENDAT/yeast
    selection=method
    out=primseg_$primseg
    cmd="$GENBRO/src/plotFeature.R -i $genome --coor $coor -r $range  -s $selection -S $settings  -f png -v -o $gene --height 4 --width 6"
    echo $cmd
    $cmd 
fi