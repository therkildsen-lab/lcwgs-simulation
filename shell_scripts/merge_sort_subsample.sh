#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/sim/rep_'$REP_ID'/'
for k in {1..2}; do
  l=$(($k+100))
  ## merge
  samtools merge $OUTDIR'sample_'$k'.bam' \
  $OUTDIR'derived_'$k'.bam' \
  $OUTDIR'derived_'$l'.bam'
  ## sort
  samtools sort -o $OUTDIR'sample_'$k'_sorted.bam' $OUTDIR'sample_'$k'.bam'
  ## subsample
  samtools view \
  -s 0.05 \
  -b $OUTDIR'sample_'$k'_sorted.bam' \
  > sample_'$k'_sorted_1x.bam
  ## delete intermediate files
  #rm $OUTDIR'derived_'$k'.bam' 
  #rm $OUTDIR'derived_'$j'.bam' 
  rm $OUTDIR'sample_'$k'.bam' 
done
