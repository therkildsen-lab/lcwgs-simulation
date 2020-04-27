#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=18
## Generate sam files
COUNT=0
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      /workdir/programs/art_bin_MountRainier/art_illumina \
      -ss HS25 \
      -sam \
      -i $OUT_DIR'fasta/p'$j'_sample'$i'_genome'$k'.fasta' \
      -p \
      -na \
      -l 150 \
      -f 10 \
      -m 500 \
      -s 75 \
      -rs $(($REP_ID*10000+$i)) \
      -o $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
        wait
        COUNT=0
      fi
    done
  done
done
wait
## Generate bam files
COUNT=0
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      samtools view -bS -F 4 $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.sam' \
      > $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.bam' &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
        wait
        COUNT=0
      fi
    done
  done
done
wait
## Sort bam files
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      rm $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.sam'
      mv $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'1.fq' $OUT_DIR'fastq/'
      mv $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'2.fq' $OUT_DIR'fastq/'
    done
  done
done
