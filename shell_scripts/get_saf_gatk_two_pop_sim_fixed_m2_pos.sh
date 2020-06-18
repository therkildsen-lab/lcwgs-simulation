#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=16
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
      for POP in {1,2}; do
      ## Get SAF
      /workdir/programs/angsd0.931/angsd/angsd \
      -b $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -out $BASE_DIR'angsd_gatk/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -dosaf 1 -GL 2 -doGlf 2 -doMajorMinor 5 \
      -doCounts 1 -doDepth 1 -dumpCounts 1 \
      -P 2 \
      -setMinDepth `awk "BEGIN {print $SAMPLE_SIZE*$COVERAGE}"` \
      -minQ 20 & 
      ## Submit two jobs at a time
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_JOB_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done
