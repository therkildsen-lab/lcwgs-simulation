#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=2
COUNT=0
# Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)				
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
		/workdir/programs/angsd0.931/angsd/misc/realSFS \
    $BASE_DIR$SUBDIR'/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    $BASE_DIR$SUBDIR'/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    -P 8 \
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' & 
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
# Get alpha_beta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
		/workdir/programs/angsd0.931/angsd/misc/realSFS fst index \
    $BASE_DIR$SUBDIR'/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    $BASE_DIR$SUBDIR'/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    -P 8 \
    -sfs $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' \
    -fstout $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta' &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
# Get alpha_beta.txt
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
		/workdir/programs/angsd0.931/angsd/misc/realSFS fst print \
    $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.fst.idx' \
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' &
  done
done
wait
# Get fst
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
		awk '{ print $0 "	" $3 / $4 }' \
    $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' \
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.fst' &
  done
done
