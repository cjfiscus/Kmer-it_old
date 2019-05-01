#!/bin/bash
# kmer-it multiple run aggregator
# cjfiscus
# 2019-04-04

# read in params
source ./params

# id samples with multiple seq runs
LST=$(cut -f2 "$SEQ_LIST" | tail -n+2 | sort | uniq -d)

cd "$OUT_DIR"

# for each sample
for i in $LST
do
	# get list of files to cat
	FILES=$(awk -v i="$i" '$2 == i' "$SEQ_LIST" | cut -f1)

	echo "$FILES"

	# cat files
	zcat $(sed 's/$/.txt.gz/g' <<< "$FILES") | sort -k1 > temp.txt
	
	# aggregate counts
	awk -F"\t" '{a[$1]+=$2;}END{for(i in a)print i"\t"a[i];}' temp.txt | gzip > "$OUT_DIR"/"$i".txt.gz 

	# cleanup
	rm temp.txt
	rm $(sed 's/$/.txt.gz/g' <<< "$FILES")
done 

