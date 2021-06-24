#!/bin/bash

## This script is used to quality filter and trim poly g tails. It can process both paired end and single end data. 
BAMLIST=$1 # Path to a list of merged, deduplicated, and overlap clipped bam files. Full paths should be included. An example of such a bam list is /workdir/cod/greenland-cod/sample_lists/bam_list_1.tsv
BASEDIR=$2 # Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
REFERENCE=$3 # Path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
REFNAME=$4 # Reference name to add to output files, e.g. gadMor2

## Loop over each sample
#for SAMPLEBAM in `cat $BAMLIST`; do

#if [ -e $SAMPLEBAM'.bai' ]; then
#	echo "the file already exists"
#else
	## Index bam files
#	samtools index $SAMPLEBAM
#fi
#
#done

## Realign around in-dels
# This is done across all samples at once

## Create list of potential in-dels
if [ ! -f $BASEDIR'bam/all_samples_for_indel_realigner.intervals' ]; then
	java -Xmx40g -jar /programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
	   -T RealignerTargetCreator \
	   -R $REFERENCE \
	   -I $BAMLIST \
	   -o $BASEDIR'bam/all_samples_for_indel_realigner.intervals' \
	   -drf BadMate
fi

## Run the indel realigner tool
java -Xmx40g -jar /programs/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R $REFERENCE \
   -I $BAMLIST \
   -targetIntervals $BASEDIR'bam/all_samples_for_indel_realigner.intervals' \
   --consensusDeterminationModel USE_READS  \
   --nWayOut _realigned.bam
