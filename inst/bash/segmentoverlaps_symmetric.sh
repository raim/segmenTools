#!/bin/bash

## calculate overlap statistics between two segmentations
## called from the segmentJaccard_bed function in segmenTools,
## but usable as a standalone script.

query=$1
target=$2
gidx=$3 # genome index file
PERM=$4

#query=query.bed
#target=target.bed
#gidx=genome.idx
#PERM=5

## get segment classes
Qtypes=($(cut -f 5 $query | sort | uniq))
qlen=${#Qtypes[@]}

 
>&2 echo GENERATING $PERM PERMUTATIONS OF QUERY $query $qlen "${Qtypes[*]}"

## NOTE: re-using regenerating permutations
## TODO: make usage saver

## generate simple genome order list
if [ ! -f genome.txt ]; then
    cut -f 1 $gidx > genome.txt
fi

pfile=`basename $query | sed 's/.bed//g'`
mkdir -p tmp
pfile=tmp/$pfile

## generate randomized queries
start=1
for (( i=$start; i<=$PERM; i++ )); do
    rfile=${pfile}_random_${i}.bed
    tfile=${pfile}_random.bed
    if [ ! -f "$rfile" ]; then
	grep -P "\t\\+$" $query | bedtools shuffle -i - -g $gidx -seed $i -noOverlapping -allowBeyondChromEnd > $tfile
	grep -P "\t\\-$" $query | bedtools shuffle -i - -g $gidx -seed $i -noOverlapping -allowBeyondChromEnd  >> $tfile
	bedtools sort -i $tfile -faidx $gidx > $rfile
    else
	>&2 echo $rfile exists
    fi
done

>&2 echo "COUNTING $qlen TOTALS"


## print header file
echo -e "query\ttarget\tintersect\tunion\tjaccard\tcount\trandom\tpermutations"

## counting totals with awk
for ((i=0; i<qlen; i++)); do
    Q=${Qtypes[$i]}
    >&2 echo "counting queries"
    pQ="\t${Q}\t"
    cnt=`grep -P $pQ $query | wc -l`
    ## calculate union via bedtools merge 
    len=`grep -P $pQ $query | bedtools sort -i - | bedtools merge -s -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "$Q\t\t\t$len\t\t$cnt\t\t"
done
for ((i=0; i<qlen; i++)); do
    T=${Qtypes[$i]}
    >&2 echo "counting targets $T"
    pT="\t${T}\t"
    cnt=`grep -P $pT $target | wc -l`
    ## calculate union via bedtools merge 
    len=`grep -P $pT $target | bedtools sort -i - | bedtools merge -s -i - | awk '{print $3-$2}' - | awk -F',' '{sum+=$0;} END{print sum;}' -`
    echo -e "\t$T\t\t$len\t\t$cnt\t\t"
done

## loop through segment types
## TODO: half loop, i=1:n and j=i:n,
##       and run in both directions, sum U and I, and calculate common J
#for i in "${array[@]}"
for ((i=0; i<qlen; i++)); do
    Q=${Qtypes[$i]}
    #echo $i ${Qtypes[$i]}
    for ((j=i; j<qlen; j++)); do
	T=${Qtypes[$j]}

	>&2 echo QUERY CLASS $i $Q v TARGET CLASS $j $T

	## grep pattern
	pQ="\t${Q}\t"
	## use randomized target name for temporary files
	qbed=${target}_Q.bed
	## grep current query
	grep -P $pQ $query  > $qbed
	
	## as above, but for target
	pT="\t${T}\t"
	tbed=${target}_T.bed
	grep -P $pT $target > $tbed
	
	## jaccard index
	jaccard=`bedtools jaccard -a $qbed -b $tbed -s -nonamecheck| grep -v intersection`
	I1=`echo $jaccard | cut -f 1 -d " "`
	U1=`echo $jaccard | cut -f 2 -d " "`
	## count of overlapping
	count1=`bedtools intersect -a $qbed -b $tbed -s  -nonamecheck|wc -l`

	## REVERSE
	grep -P $pT $query  > $qbed
	grep -P $pQ $target > $tbed

	## jaccard index
	jaccard=`bedtools jaccard -a $qbed -b $tbed -s -nonamecheck| grep -v intersection`
	I2=`echo $jaccard | cut -f 1 -d " "`
	U2=`echo $jaccard | cut -f 2 -d " "`
	## count of overlapping
	count2=`bedtools intersect -a $qbed -b $tbed -s  -nonamecheck|wc -l`

	## sum up I/U/count values and calculate J with awk
	I=$(awk -v i1="$I1" -v i2="$I2" 'BEGIN{print i1+i2}')
	U=$(awk -v u1="$U1" -v u2="$U2" 'BEGIN{print u1+u2}')
	J=0
	isntzero=$(echo "$U != 0" | bc -l)
	if ((isntzero == 1)); then
	    J=$(awk -v i="$I" -v u="$U" 'BEGIN{print i/u}')
	fi
	count=$(awk -v c1="$count1" -v c2="$count2" 'BEGIN{print c1+c2}')

	## permutation p-value
	cnt=0
	tot=0
	start=1
	for (( p=$start; p<=$PERM; p++ )) 
	do
	    let tot++;
	    rfile=${pfile}_random_${p}.bed

	    grep -P $pT $target  > $tbed
	    jaccard=`grep -P $pQ $rfile | bedtools jaccard -a - -b $tbed -s -nonamecheck | grep -v intersection`
	    Ir1=`echo $jaccard | cut -f 1 -d " "`
	    Ur1=`echo $jaccard | cut -f 2 -d " "`

	    ## REVERSE
	    grep -P $pQ $target  > $tbed
	    jaccard=`grep -P $pT $rfile | bedtools jaccard -a - -b $tbed -s -nonamecheck | grep -v intersection`
	    Ir1=`echo $jaccard | cut -f 1 -d " "`
	    Ur1=`echo $jaccard | cut -f 2 -d " "`

	    ## sum up I/U/count values and calculate J with awk
	    Ir=$(awk -v i1="$Ir1" -v i2="$Ir2" 'BEGIN{print i1+i2}')
	    Ur=$(awk -v u1="$Ur1" -v u2="$Ur2" 'BEGIN{print u1+u2}')
	    Jr=0
	    isntzero=$(echo "$Ur != 0" | bc -l)
	    if ((isntzero == 1)); then
		Jr=$(awk -v i="$Ir" -v u="$Ur" 'BEGIN{print i/u}')
	    fi
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
