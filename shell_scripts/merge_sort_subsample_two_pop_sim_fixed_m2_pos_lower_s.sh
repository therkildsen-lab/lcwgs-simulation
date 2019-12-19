#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_'$REP_ID'/'
N_CORE_MAX=40
## merge
COUNT=0
for k in {1..160}; do
  for i in {1..2}; do
    samtools merge $OUT_DIR'bam/p'$i'_sample_'$k'.bam' \
    $OUT_DIR'bam/p'$i'_derived_'$k'_1.bam' \
    $OUT_DIR'bam/p'$i'_derived_'$k'_2.bam' &
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
for k in {1..160}; do
  for i in {1..2}; do
    samtools sort -o $OUT_DIR'bam/p'$i'_sample_'$k'_sorted.bam' $OUT_DIR'bam/p'$i'_sample_'$k'.bam' &
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
for j in {0.25,0.5,1,2,4,8}; do
  for k in {1..160}; do
    for i in {1..2}; do
      samtools view \
      -s `awk -v j=$j 'BEGIN { print j / 20 }'` \
      -b $OUT_DIR'bam/p'$i'_sample_'$k'_sorted.bam' \
      > $OUT_DIR'bam/p'$i'_sample_'$k'_sorted_'$j'x.bam' &
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
for k in {1..160}; do
  for i in {1..2}; do
    l=$(($k+160))
    #rm $OUT_DIR'bam/p'$i'_derived_'$k'.bam' 
    #rm $OUT_DIR'bam/p'$i'_derived_'$l'.bam' 
    rm $OUT_DIR'bam/p'$i'_sample_'$k'.bam'
  done
done
