#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=42
## merge
COUNT=0
for k in {1..5}; do
  for i in {1..9}; do
    samtools merge $OUT_DIR'bam/p'$i'_sample'$k'.bam' \
    $OUT_DIR'bam/p'$i'_sample'$k'_genome1.bam' \
    $OUT_DIR'bam/p'$i'_sample'$k'_genome2.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## sort
COUNT=0
for k in {1..5}; do
  for i in {1..9}; do
    samtools sort -o $OUT_DIR'bam/p'$i'_sample'$k'_sorted.bam' $OUT_DIR'bam/p'$i'_sample'$k'.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## subsample
COUNT=0
for j in {0.125,0.25,0.5,1,2,4}; do
  for k in {1..5}; do
    for i in {1..9}; do
      samtools view \
      -s `awk -v j=$j 'BEGIN { print j / 20 }'` \
      -b $OUT_DIR'bam/p'$i'_sample'$k'_sorted.bam' \
      > $OUT_DIR'bam/p'$i'_sample'$k'_sorted_'$j'x.bam' &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done
wait
## delete intermediate files
for k in {1..5}; do
  for i in {1..9}; do
    #rm $OUT_DIR'bam/p'$i'_sample'$k'_genome1.bam' 
    #rm $OUT_DIR'bam/p'$i'_sample'$k'_genome2.bam' 
    rm $OUT_DIR'bam/p'$i'_sample'$k'.bam'
  done
done
