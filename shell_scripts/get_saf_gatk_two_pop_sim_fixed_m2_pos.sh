#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=4
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
      for POP in {1,2}; do
      ## Get SAF
      /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR'angsd_gatk/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -doSaf 1 \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL 2 \
      -P 4 \
      -doCounts 1 \
      -setMinDepth `awk "BEGIN {print $SAMPLE_SIZE*$COVERAGE}"` & 
      ## Submit four jobs at a time
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_JOB_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done
