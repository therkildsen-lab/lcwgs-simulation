#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
for SAMPLE_SIZE in {5,10,20,40,80}; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
      ## Run PCAngsd
      python /workdir/programs/pcangsd/pcangsd.py \
      -beagle $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.beagle.gz' \
      -minMaf 0.05 \
      -threads 16 \
      -o $BASE_DIR'angsd/pcagnsd_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x'
  done
done
