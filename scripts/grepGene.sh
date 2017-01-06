#!/usr/bin/env bash

## TODO: parse arguments

## gene: standard gene name, IME4, IME1
gene=$1 
features=/data/yeast/feature_R64-1-1_20110208_withclusters.csv
sgtype=$YEASTSEQ/segmentation/segmentTest/20161228/annotation/T.ash_D.dcash_K.16_S.icor_E.3_M.175_

seg=${sgtype}genes.csv
sas=${sgtype}antisense.csv
id=`sed 's/"//g' $features|grep gene | grep -P "\t${gene}\t"| cut -f 2`
sid=`grep $id $seg | cut -f 1|sed 's/.*\://'`
primseg=`echo $sid|sed s/_.*//`
asid=`grep -P "^$sid\t" $sas| cut -f 2`
asprimseg=`echo $asid|sed s/_.*//`


echo GENE $gene 
echo PRIMSEG $primseg SEGMENT $sid 
echo ANTISENSE $asprimseg SEGMENT $asid 
