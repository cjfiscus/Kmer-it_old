#!/bin/bash -l

# Kmer-it
# cjfiscus
# 2019-03-25

##### PIPELINE #####
# read in arguments from params file
source $1

# determine sample name
NAME=$(head -n "$2" "$SEQ_LIST" | tail -n 1 | cut -f1)
echo "$NAME"

# determine if seq run is SE or PE
FILE=$(head -n "$2" "$SEQ_LIST" | tail -n 1 | cut -f3)
echo "$FILE"

if [[ "$FILE" == *";"* ]] ; then
        echo "PE library detected"
	LIBTYPE="PE" # paired end
else
	echo "SE library detected"
        LIBTYPE="SE" # single end
fi

# work in temp directory
TEMP_DIR="$TEMP_DIR"/"$NAME"
mkdir -pv "$TEMP_DIR"
cd "$TEMP_DIR"

# make outdir for mapstats
mkdir "$OUT_DIR"/mapstats

if [ $LIBTYPE == "PE" ]
then # paired end 
        # Download files
	INDEX=1 # 1 is forward, 2 is reverse
	for i in $(echo $FILE | tr ";" "\n")
	do	
		echo "downloading" "$i" 
		wget "$i" -O "$NAME"_"$INDEX".fastq.gz
		INDEX=$((INDEX + 1))
	done

	# check MD5sums
	INDEX=1
	MD5SUMS=$(head -n "$2" "$SEQ_LIST" | tail -n 1 | cut -f4)
	echo "$MD5SUMS"
	if [ -z "$MD5SUMS" ]
	then 
		echo "skipping MD5sum check..."
	else
		echo "checking MD5sums..."

		for i in $(echo "$MD5SUMS" | tr ";" "\n")
		do 
			echo "$i" "$NAME"_"$INDEX".fastq.gz >> chk.md5
    			INDEX=$((INDEX + 1))
		done	

		if md5sum --status -c chk.md5; then
			# continue		
			echo "SUMS OK"
		else
			# stop script
			echo "SUMS BAD"
			exit 1
		fi
	fi

	if [[ $RUN_TRIM = "yes" ]]
	then 
		# Quality/Adapter trimming
		echo "trimming with trimmomatic..."
		java -jar $TRIMMOMATIC PE -threads "$THREADS" \
		"$NAME"_1.fastq.gz "$NAME"_2.fastq.gz \
		"$NAME"_1_trimmed_paired.fq.gz "$NAME"_1_unpaired.fq.gz \
		"$NAME"_2_trimmed_paired.fq.gz "$NAME"_2_unpaired.fq.gz \
		ILLUMINACLIP:"$ADAPTERSPE":2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

	else
		echo "skipping trimming..."
		mv "$NAME"_1.fastq.gz "$NAME"_1_trimmed_paired.fq.gz
		mv "$NAME"_2.fastq.gz "$NAME"_2_trimmed_paired.fq.gz
	fi

	if [ -z "$O_GENOME" ]
	then 
		echo "No mapping to organellar genome"
		zcat "$NAME"_*_trimmed_paired.fq.gz > "$NAME".unmapped.fq
	
	else
		# map to organellar genome
		echo "mapping to organellar genome with bwa..."
		bwa mem -t "$THREADS" -M $O_GENOME "$NAME"_1_trimmed_paired.fq.gz \
			"$NAME"_2_trimmed_paired.fq.gz > "$NAME"_org.sam
	fi

else # single end 
	# Download file
	echo "downloading" "$FILE"	
	wget "$FILE" -O "$NAME".fastq.gz
	
	# check MD5sum
	MD5SUMS=$(head -n "$2" $SEQ_LIST | tail -n 1 | cut -f4)
	if [ -z "$MD5SUMS" ]
	then 
		echo "skipping MD5sum check..."
	else
		echo "verifying checksums..."
		echo "$MD5SUMS" "$NAME".fastq.gz >> chk.md5
	
		if md5sum --status -c chk.md5; then
        		# continue
        		echo "SUMS OK"
        	else 
        		# stop script
        		echo "SUMS BAD"
        		exit 1
        	fi
	fi 

	if [[ $RUN_TRIM = "yes" ]]
	then 
		# Quality/Adapter trimming
		echo "trimming with trimmomatic..."
		java -jar $TRIMMOMATIC SE -threads $THREADS \
		"$NAME".fastq.gz "$NAME"_trimmed.fq.gz \
		ILLUMINACLIP:"$ADAPTERSSE":2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

	else
		echo "skipping trimming..."
		mv "$NAME".fastq.gz "$NAME"_trimmed.fq.gz
	fi

	if [ -z "$O_GENOME" ]
	then
		echo "No mapping to organellar genome" 
		zcat "$NAME"_trimmed.fq.gz > "$NAME".unmapped.fq
	else 
		# map to organellar genome
		echo "mapping to organellar genome with bwa..."
		bwa mem -t "$THREADS" -M $O_GENOME "$NAME"_trimmed.fq.gz  > "$NAME"_org.sam

	fi 
fi

# sam to sorted bam
if [ -n "$O_GENOME" ]
then 
samtools view -bS "$NAME"_org.sam | samtools sort -T temp_Pt - -o "$NAME"_org.bam
samtools flagstat "$NAME"_org.bam > $OUT_DIR/mapstats/"$NAME"_org_mapstats.txt

echo "extracting unmapped reads..."
samtools view -f4 -b "$NAME"_org.bam > "$NAME".unmapped.bam

# export these unmapped reads
bedtools bamtofastq -i "$NAME".unmapped.bam -fq "$NAME".unmapped.fq
fi

# Count K-mers in reads that did not map to organelles 
echo "counting K-mers with jellyfish"
jellyfish count -C -m "$K" -s 500M -t "$THREADS" -o "$NAME".jf "$NAME".unmapped.fq 
jellyfish dump -tc "$NAME".jf | gzip > $OUT_DIR/"$NAME".txt.gz


# Delete temp folder (if enabled)
if [[ $RM_TEMP_DIR = "yes" ]]
then 
	rm -r "$TEMP_DIR"
fi
