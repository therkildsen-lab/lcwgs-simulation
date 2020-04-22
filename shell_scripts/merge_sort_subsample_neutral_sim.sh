#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/neutral_sim/}

OUT_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=30
## merge
COUNT=0
for k in {1..1000}; do
  i=$((2*$k-1))
  j=$((2*$k))
  ## merge
  samtools merge $OUT_DIR'bam/sample_'$k'.bam' \
  $OUT_DIR'bam/derived_'$i'.bam' \
  $OUT_DIR'bam/derived_'$j'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## sort
COUNT=0
for k in {1..1000}; do
  samtools sort -o $OUT_DIR'bam/sample_'$k'_sorted.bam' $OUT_DIR'bam/sample_'$k'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## subsample
COUNT=0
for k in {1..1000}; do
  for j in {0.25,0.5,1,2,4,8}; do
    samtools view \
    -s `awk -v j=$j 'BEGIN { print j / 20 }'` \
    -b $OUT_DIR'bam/sample_'$k'_sorted.bam' \
    > $OUT_DIR'bam/sample_'$k'_sorted_'$j'x.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
for k in {1..1000}; do
  i=$((2*$k-1))
  j=$((2*$k))
  ## delete intermediate files
  #rm $OUT_DIR'bam/derived_'$i'.bam' 
  #rm $OUT_DIR'bam/derived_'$j'.bam' 
  rm $OUT_DIR'bam/sample_'$k'.bam' 
done
