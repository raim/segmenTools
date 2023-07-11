#!/bin/bash

## calculate overlap statistics between two segmentations
## called from the segmentJaccard_bed function in segmenTools,
## but usable as a standalone script.

query=$1
target=$2
gidx=$3
PERM=$4
keep=$5

## awk function for division
calc(){ awk "BEGIN { print "$*" }"; }



## generate simple genome order list
cut -f 1 $gidx > genome.txt

## get segment classes
Qtypes=`cut -f 5 $query | sort | uniq`
Ttypes=`cut -f 5 $target | sort | uniq`


>&2 echo GENERATING $PERM PERMUTATIONS OF QUERY $query $Qtypes

## NOTE: storing permutations in main dir during usage,
## to save time by not regenerating permutations
## DELETE once intended for broader use.

##TMPDIR=$(dirname $(mktemp -u))
pfile=`basename $query | sed 's/.bed//g'`

start=1
for (( i=$start; i<=$PERM; i++ )); do
    rfile=${pfile}_random_${i}.bed
    if [ ! -f "$rfile" ]; then
	bedtools shuffle -i $query -g $gidx -seed $i  | bedtools sort -i - -faidx  genome.txt > $rfile
    else
	>&2 echo $rfile exists
    fi
done


>&2 echo CALCULATING OVERLAPS TO TARGET $target $Ttypes

## print header file
echo -e "query\ttarget\tintersect\tunion\tjaccard\tcount\trandom\tpermutations"

## counting totals with awk
for Q in $Qtypes; do
    pQ="\t${Q}\t"
    cnt=`grep -P $pQ $query | wc -l`
    ## calculate union via bedtools merge (
    len=`grep -P $pQ $query | bedtools merge -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "$Q\t\t\t$len\t\t$cnt\t\t"
done
for T in $Ttypes; do
    pT="\t${T}\t"
    cnt=`grep -P $pT $target | wc -l`
    ## calculate union via bedtools merge (
    len=`grep -P $pT $target | bedtools merge -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "\t$T\t\t$len\t\t$cnt\t\t"
done

## loop through segment types
for Q in $Qtypes; do
    for T in $Ttypes; do

	>&2 echo QUERY CLASS $Q v TARGET CLASS $T

	## grep pattern
	pQ="\t${Q}\t"
	pT="\t${T}\t"
	
	grep -P $pQ $query  > Q.bed
	grep -P $pT $target > T.bed
	
	## jaccard index
	jaccard=`bedtools jaccard -a Q.bed -b T.bed -s | grep -v intersection`
	I=`echo $jaccard | cut -f 1 -d " "`
	U=`echo $jaccard | cut -f 2 -d " "`
	J=`echo $jaccard | cut -f 3 -d " "`
	
	## count of overlapping
	count=`bedtools intersect -a Q.bed -b T.bed -s |wc -l`
		
	## permutation p-value
	cnt=0
	tot=0
	start=1
	for (( i=$start; i<=$PERM; i++ )) 
	do
	    let tot++;
	    rfile=${pfile}_random_${i}.bed
	    Jr=`grep -P $pQ $rfile | bedtools jaccard -a - -b T.bed -s | grep -v intersection | cut -f 3`
	    if [ $Jr \> $J ]; then let cnt++; fi
	done
	## avoid calculation in bash, unless sure of it.
	##pvalue=`calc $cnt/$tot`

	## print result
	echo -e "$Q\t$T\t$I\t$U\t$J\t$count\t$cnt\t$tot" 
    done
done 
## partially clean up; note that main cleaning happens via the
## segmentJaccard_bed function in segmenTools
\rm -f genome.txt