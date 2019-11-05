#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
for k in {1..160}; do
  l=$(($k+160))
  ## merge
  samtools merge $OUT_DIR'bam/sample_'$k'.bam' \
  $OUT_DIR'bam/derived_'$k'.bam' \
  $OUT_DIR'bam/derived_'$l'.bam'
  ## sort
  samtools sort -o $OUT_DIR'bam/sample_'$k'_sorted.bam' $OUT_DIR'bam/sample_'$k'.bam'
  ## subsample
  for j in {0.25,0.5,1,2,4,8}; do
    samtools view \
    -s `awk -v j=$j 'BEGIN { print j / 20 }'` \
    -b $OUT_DIR'bam/sample_'$k'_sorted.bam' \
    > $OUT_DIR'bam/sample_'$k'_sorted_'$j'x.bam'
  done
  ## delete intermediate files
  #rm $OUT_DIR'bam/derived_'$k'.bam' 
  #rm $OUT_DIR'bam/derived_'$l'.bam' 
  rm $OUT_DIR'bam/sample_'$k'.bam' 
done
