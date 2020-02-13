#!/bin/bash
REP_ID=$1
IN_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
OUT_DIR='/workdir/lcwgs-simulation/neutral_sim_uneven_input/rep_'$REP_ID'/'
N_CORE_MAX=28

## subsample
COUNT=0
for i in {5,10,20,40,80,160}; do
  for k in $(seq 1 $i); do
    ROW_INDEX=$k
    COL_INDEX=0
    for j in {0.25,0.5,1,2,4,8}; do
      COL_INDEX=$(( COL_INDEX + 1 ))
      COVERAGE=`cat $OUT_DIR'misc/coverage_matrix_'$i'.tsv' | head -n $ROW_INDEX | tail -n 1 | cut -f $COL_INDEX`
      echo $COVERAGE
      samtools view \
      -s `awk -v j=$COVERAGE 'BEGIN { print j / 20 }'` \
      -b $IN_DIR'bam/sample_'$k'_sorted.bam' \
      > $OUT_DIR'bam/sample_'$k'_sorted_'$i'_'$j'x.bam' &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
        wait
        COUNT=0
      fi
    done
  done
done
