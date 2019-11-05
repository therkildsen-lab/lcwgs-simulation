#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
for i in {1..320}; do
  /workdir/programs/art_bin_MountRainier/art_illumina \
  -ss HS25 \
  -sam \
  -i $OUT_DIR'fasta/derived_'$i'.fasta' \
  -p \
  -na \
  -l 150 \
  -f 10 \
  -m 500 \
  -s 75 \
  -rs $(($REP_ID*10000+$i)) \
  -o $OUT_DIR'bam/derived_'$i
samtools view -bS -F 4 $OUT_DIR'bam/derived_'$i'.sam' > $OUT_DIR'bam/derived_'$i'.bam'
rm $OUT_DIR'bam/derived_'$i'.sam'
mv $OUT_DIR'bam/derived_'$i'1.fq' $OUT_DIR'fastq/'
mv $OUT_DIR'bam/derived_'$i'2.fq' $OUT_DIR'fastq/'
done
