#!/bin/bash

## calculate overlap statistics between two segmentations
## called from the segmentJaccard_bed function in segmenTools,
## but usable as a standalone script.

query=$1
target=$2  
gidx=$3    # genome index file
PERM=$4    # number of permutations
overlap=$5 # allow overlaps? pass yes 

## get segment classes
Qtypes=`cut -f 5 $query | sort | uniq`
Ttypes=`cut -f 5 $target | sort | uniq`


>&2 echo GENERATING $PERM PERMUTATIONS OF QUERY $query $Qtypes
if [ "$overlap"  = "yes" ]; then
    >&2 echo -e "\tALLOWING OVERLAPS IN PERMUTATIONS"
fi

## NOTE: re-using regenerating permutations
## TODO: make usage saver

## generate simple genome order list
if [ ! -f genome.txt ]; then
    cut -f 1 $gidx > genome.txt
fi

pfile=`basename $query | sed 's/.bed//g'`
mkdir -p tmp
pfile=tmp/$pfile

## TODO: in bedtools v2.27.1 shuffle doesn't finish with -noOverlapping option
## (>10k input segments), but it works with -allowBeyondChromEnd.
## HOWEVER, this will also no allow overlapping features from opposite strands,
## NOTE: e.g. eccDNA_all : all segment clusters become significant with p=0
## when randomization occurs with -noOverlapping -allowBeyondChromEnd
## NOTE: shuffling each strand separately with -noOverlapping option, to allow
## antisense overlaps; most closely reflects original data and internal
## randomization, but also requires the -allowBeyondChromEnd option to finish.

## generate randomized queries
start=1
for (( i=$start; i<=$PERM; i++ )); do
    rfile=${pfile}_random_${i}.bed
    tfile=${pfile}_random.bed
    if [ ! -f "$rfile" ]; then
	if [ "$overlap"  = "yes" ]; then
	    ## allowing overlaps:
	    bedtools shuffle -i $query -g $gidx -seed $i | bedtools sort -i - -faidx $gidx > $rfile
	else
	    ## no overlaps, by strand!
	    grep -P "\t\\+$" $query | bedtools shuffle -i - -g $gidx -seed $i -noOverlapping -allowBeyondChromEnd > $tfile
	    grep -P "\t\\-$" $query | bedtools shuffle -i - -g $gidx -seed $i -noOverlapping -allowBeyondChromEnd  >> $tfile
	    bedtools sort -i $tfile -faidx $gidx > $rfile
	fi
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
    ## calculate union via bedtools merge 
    len=`grep -P $pQ $query | bedtools sort -i - | bedtools merge -s -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "$Q\t\t\t$len\t\t$cnt\t\t"
done
for T in $Ttypes; do
    pT="\t${T}\t"
    cnt=`grep -P $pT $target | wc -l`
    ## calculate union via bedtools merge 
    len=`grep -P $pT $target | bedtools sort -i - | bedtools merge -s -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "\t$T\t\t$len\t\t$cnt\t\t"
done

## loop through segment types
for Q in $Qtypes; do

    ## grep pattern
    pQ="\t${Q}\t"
    ## use randomized target name for temporary files
    qbed=${target}_Q.bed
    ## grep current query
    grep -P $pQ $query  > $qbed

    for T in $Ttypes; do

	>&2 echo QUERY CLASS $Q v TARGET CLASS $T

	## as above, but for target
	pT="\t${T}\t"
	tbed=${target}_T.bed
	grep -P $pT $target > $tbed
	
	## jaccard index
	jaccard=`bedtools jaccard -a $qbed -b $tbed -s -nonamecheck| grep -v intersection`
	I=`echo $jaccard | cut -f 1 -d " "`
	U=`echo $jaccard | cut -f 2 -d " "`
	J=`echo $jaccard | cut -f 3 -d " "`

	## count of overlapping
	count=`bedtools intersect -a $qbed -b $tbed -s  -nonamecheck|wc -l`
		
	## permutation p-value
	cnt=0
	tot=0
	start=1
	for (( i=$start; i<=$PERM; i++ )) 
	do
	    let tot++;
	    rfile=${pfile}_random_${i}.bed
	    Jr=`grep -P $pQ $rfile | bedtools jaccard -a - -b $tbed -s -nonamecheck | grep -v intersection | cut -f 3`
	    ## is Jr >= J ?
	    r=$(awk -v j="$J" -v r="$Jr" 'BEGIN{print (j<=r)?1:0}')
	    if [ 1 -eq $r ]; then
		let cnt++;
	    fi
	done

	## print result
	echo -e "$Q\t$T\t$I\t$U\t$J\t$count\t$cnt\t$tot" 
    done
done 
