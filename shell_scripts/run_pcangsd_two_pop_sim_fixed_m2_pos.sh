#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=3
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    ## Run PCA
    python /workdir/programs/pcangsd/pcangsd.py \
    -beagle $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.beagle.gz' \
    -selection \
    -sites_save \
    -minMaf 0.05 \
    -threads 10 \
    -iter 200 \
    -maf_iter 200 \
    -o $BASE_DIR'angsd/pcagnsd_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' &
    ## Submit three jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
