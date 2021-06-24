#!/bin/bash

BAMLIST=$1 # List of base file names for the sequenced samples (without extensions _R1.fastq.gz and _R2.fastq.gz), e.g. /workdir/Cod/Greenland/SampleLists/SampleList_PE.txt
BASEDIR=$2
NREADS=$3
SAMPLEINFO=$4
COVERAGE=$5

for SAMPLEFILE in `cat $BAMLIST`; do

SAMPLE=`grep -P "${SAMPLEFILE}\t" $SAMPLEINFO | cut -f 4`

frac=$(samtools idxstats $BASEDIR$SAMPLEFILE | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=$NREADS/total; if (frac > 1) {print 1} else {print frac}}')

samtools view -bs $frac $SAMPLEFILE > $BASEDIR$SAMPLE'_'$COVERAGE'.bam'

done