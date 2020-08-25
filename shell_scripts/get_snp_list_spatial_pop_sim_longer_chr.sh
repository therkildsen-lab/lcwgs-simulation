#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=42
COUNT=0
## Create SNP lists
for SAMPLE_SIZE in 5; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
    gunzip -c $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > \
    $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
    ## Submit a certain number of jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
COUNT=0
## Index the SNP lists
for SAMPLE_SIZE in 5; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
    /workdir/programs/angsd0.931/angsd/angsd sites index \
    $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
    ## Submit a certain number of jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
