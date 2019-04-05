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
mkdir "$TEMP_DIR"
cd "$TEMP_DIR"

if [ $LIBTYPE == "PE" ]
then # paired end 
        # Download files
	INDEX=1 # 1 is forward, 2 is reverse
	for i in $(echo $FILE | tr ";" "\n")
	do	
		echo "downloading" "$i" 
		axel -n "$THREADS" "$i" -o "$NAME"_"$INDEX".fastq.gz
		#wget "$i" -O "$NAME"_"$INDEX".fastq.gz
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

	if [ -z "$REF_GENOME" ]
	then 
		echo "No mapping to reference genome"

	else
		# map to reference genome
		echo "mapping to genome with bwa..."
		bwa mem -t "$THREADS" -M $REF_GENOME "$NAME"_1_trimmed_paired.fq.gz \
			"$NAME"_2_trimmed_paired.fq.gz > "$NAME"_gen.sam

		if [[ $REP_ASSEM = "yes"]]
		then 
			# calculate insert size metrics
        		java -jar $PICARD CollectInsertSizeMetrics \
        		I="$NAME"_gen.sam \
        		O="$OUT_DIR"/"$NAME"_metrics.txt \
        		H="$NAME"_histogram.pdf \
        		M=0.5

			# parse insert size matrics
       			MEAN_INSERT_SIZE=$(head -n 8 "$OUT_DIR"/"$NAME"_metrics.txt | tail -n 1 | cut -f6)
        		STD_INSERT_SIZE=$(head -n 8 "$OUT_DIR"/"$NAME"_metrics.txt | tail -n 1 | cut -f7)

        		# prepare REPdenovo inputs
        		touch reads
        		echo "$NAME"_1.fastq.gz 1 $(echo "$MEAN_INSERT_SIZE" | awk '{print int($1+0.5)}') $(echo "$STD_INSERT_SIZE" | awk '{print int($1+0.5)}') >> reads
        		echo "$NAME"_2.fastq.gz 1 $(echo "$MEAN_INSERT_SIZE" | awk '{print int($1+0.5)}') $(echo "$STD_INSERT_SIZE" | awk '{print int($1+0.5)}') >> reads
        
        		LEN=$(head -n2 sample1.fq | tail -n 1 | wc -m | awk '{print $1 - 1}')

		fi 
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
	axel -n "$THREADS" "$FILE" -o "$NAME".fastq.gz
	#wget "$FILE" -O "$NAME".fastq.gz
	
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
		java -jar $TRIMMOMATIC SE -threads 8 \
		"$NAME".fastq.gz "$NAME"_trimmed.fq.gz \
		ILLUMINACLIP:"$ADAPTERSSE":2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

	else
		echo "skipping trimming..."
		mv "$NAME".fastq.gz "$NAME"_trimmed.fq.gz
	fi

	if [ -z "$REF_GENOME" ]
	then
		echo "No mapping to reference genome"
	else
		# map to reference genome
		echo "mapping to genome with bwa..."
		bwa mem -t "$THREADS" -M $REF_GENOME "$NAME"_trimmed.fq.gz  > "$NAME"_gen.sam
	fi 

	if [[ $REP_ASSEM = "yes" ]]
	then 
		touch reads
		echo "$NAME".fastq.gz -1 -1 -1 >> reads
		LEN=$(zless "$NAME".fastq.gz | head -n2 | tail -n 1 | wc -m | awk '{print $1 - 1}')
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
if [ -n "$REF_GENOME" ] 
then 
samtools view -bS "$NAME"_gen.sam | samtools sort -T temp_Pt - -o "$NAME"_gen.bam
samtools flagstat "$NAME"_gen.bam > $OUT_DIR/"$NAME"_gen_mapstats.txt
samtools index "$NAME"_gen.bam

# calculate coverage of ref per base and in 1kb windows
echo "calculating coverage with mosdepth..."
mosdepth -t "$THREADS" -b 1000 "$OUT_DIR"/"$NAME" "$NAME"_gen.bam
fi 

if [ -n "$O_GENOME" ]
then 
samtools view -bS "$NAME"_org.sam | samtools sort -T temp_Pt - -o "$NAME"_org.bam
samtools flagstat "$NAME"_org.bam > $OUT_DIR/"$NAME"_org_mapstats.txt

echo "extracting unmapped reads..."
samtools view -f4 -b "$NAME"_org.bam > "$NAME".unmapped.bam

# export these unmapped reads
bedtools bamtofastq -i "$NAME".unmapped.bam -fq "$NAME".unmapped.fq
fi

# Count K-mers in reads that did not map to organelles 
echo "counting K-mers with jellyfish"
jellyfish count -C -m "$K" -s 500M -t "$THREADS" -o "$NAME".jf "$NAME".unmapped.fq 
jellyfish dump -tc "$NAME".jf | gzip > $OUT_DIR/"$NAME".txt.gz

# Repeat Assembly with REPdenovo
if [[ $REP_ASSEM = "yes" ]]
then 
	echo "assembling repeats with REPdenovo"
	# make config file for REPdenovo
	touch configure
	head -n 7 "$REP_CONFIG" >> configure
	echo READ_LENGTH "$LEN" >> configure
	tail -n+9 "$REP_CONFIG" >> configure
	echo "OUTPUT_FOLDER $TEMP_DIR" >> configure
	echo "VERBOSE 1" >> configure

	# Run REPdenovo pipeline
	python "$REPDENOVO" -c All -g configure -r reads

	# save results
	gzip contigs.fa
	mv contigs.fa.gz "$OUT_DIR"/"$NAME"_contigs.fa.gz
fi 
