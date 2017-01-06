#!/usr/bin/env bash

## TODO: parse arguments

## gene: standard gene name, IME4, IME1
gene=$1 
features=/data/yeast/feature_R64-1-1_20110208_withclusters.csv
sgtype=$YEASTSEQ/segmentation/segmentTest/20161228/annotation/T.ash_D.dcash_K.16_S.icor_E.3_M.175_

seg=${sgtype}genes.csv
sas=${sgtype}antisense.csv
id=`grep $gene $features | grep gene | cut -f 2 | sed 's/"//g'`
sid=`grep $id $seg | cut -f 1`
primseg=`echo $sid|sed s/_.*//`
asid=`grep $sid $sas| cut -f 2`
asprimseg=`echo $asid|sed s/_.*//`

echo GENE $gene 
echo SEGMENT $sid PRIMSEG $primseg
echo ANTISENSE $asid PRIMSEG $ asprimseg
