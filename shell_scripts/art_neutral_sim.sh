#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/neutral_sim/}

OUT_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=30
## Generate sam files
COUNT=0
for i in {1..2000}; do
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
  -o $OUT_DIR'bam/derived_'$i &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## Generate bam files
COUNT=0
for i in {1..2000}; do
  samtools view -bS -F 4 $OUT_DIR'bam/derived_'$i'.sam' > $OUT_DIR'bam/derived_'$i'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## Sort bam files
for i in {1..2000}; do
  rm $OUT_DIR'bam/derived_'$i'.sam'
  mv $OUT_DIR'bam/derived_'$i'1.fq' $OUT_DIR'fastq/'
  mv $OUT_DIR'bam/derived_'$i'2.fq' $OUT_DIR'fastq/'
done
