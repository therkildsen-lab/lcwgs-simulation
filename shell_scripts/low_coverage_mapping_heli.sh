#!/bin/bash

## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data. 
SAMPLELIST=$1 # Path to a list of prefixes of the raw fastq files. It should be a subset of the the 1st column of the sample table. An example of such a sample list is /workdir/cod/greenland-cod/sample_lists/sample_list_pe_1.tsv
SAMPLETABLE=$2 # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table.tsv
FASTQDIR=$3 # Path to the directory where fastq file are stored. An example for the quality-filtered Greenland pe data is: /workdir/cod/greenland-cod/qual_filtered/
BASEDIR=$4 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
FASTQSUFFIX1=$5 # Suffix to fastq files. Use forward reads with paired-end data. An example for the quality-filtered Greenland paired-end data is: _adapter_clipped_qual_filtered_f_paired.fastq.gz
FASTQSUFFIX2=$6 # Suffix to fastq files. Use reverse reads with paired-end data. An example for the quality-filtered Greenland paired-end data is: _adapter_clipped_qual_filtered_r_paired.fastq.gz
MAPPINGPRESET=$7 # The pre-set option to use for mapping in bowtie2 (very-sensitive for end-to-end (global) mapping [typically used when we have a full genome reference], very-sensitive-local for partial read mapping that allows soft-clipping [typically used when mapping genomic reads to a transcriptome]
REFERENCE=$8 # Path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
REFNAME=$9 # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
for SAMPLEFILE in `cat $SAMPLELIST`; do
	
	## Extract relevant values from a table of sample, sequencing, and lane ID (here in columns 4, 3, 2, respectively) for each sequenced library
	SAMPLE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 4`
	SEQ_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 3`
	LANE_ID=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	SAMPLE_SEQ_ID=$SAMPLE_ID'_'$SEQ_ID'_'$LANE_ID
	
	## Extract data type from the sample table
	DATATYPE=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 6`
	
	## The input and output path and file prefix
	SAMPLETOMAP=$FASTQDIR$SEQ_ID
	SAMPLEBAM=$BASEDIR'bam/'$SAMPLE_SEQ_ID
	
	## Define platform unit (PU), which is the lane number
	PU=`grep -P "${SAMPLEFILE}\t" $SAMPLETABLE | cut -f 2`
	
	## Define reference base name
	REFBASENAME="${REFERENCE%.*}"
	
	## Map reads to the reference 
	# Map the paired-end reads
	if [ $DATATYPE = pe ]; then 
	# We ignore the reads that get orphaned during adapter clipping because that is typically a very small proportion of reads. If a large proportion of reads get orphaned (loose their mate so they become single-end), these can be mapped in a separate step and the resulting bam files merged with the paired-end mapped reads.
		bowtie2 -q --phred33 --$MAPPINGPRESET -p 16 -I 0 -X 1500 --fr --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -1 $SAMPLETOMAP$FASTQSUFFIX1 -2 $SAMPLETOMAP$FASTQSUFFIX2 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'
	
	# Map the single-end reads
	elif [ $DATATYPE = se ]; then
		bowtie2 -q --phred33 --$MAPPINGPRESET -p 16 --rg-id $SAMPLE_SEQ_ID --rg SM:$SAMPLE_ID --rg LB:$SAMPLE_ID --rg PU:$PU --rg PL:ILLUMINA -x $REFBASENAME -U $SAMPLETOMAP$FASTQSUFFIX1 -S $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'
	
	fi
	
	## Convert to bam file for storage
	samtools view -bS -F 4 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam' > $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam'
	rm $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.sam'
	
	## Filter the mapped reads
	# Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20) -- do we want the quality score filter??
	samtools view -h -q 20 $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'.bam' | samtools view -buS - | samtools sort -o $SAMPLEBAM'_'$DATATYPE'_bt2_'$REFNAME'_minq20_sorted.bam'
	
done
