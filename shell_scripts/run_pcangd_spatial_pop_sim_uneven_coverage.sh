#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
for SAMPLE_SIZE in {5,10,20,40,80}; do
  ## Run PCAngsd
  python2 /workdir/programs/pcangsd/pcangsd.py \
  -beagle $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_uneven_coverage.beagle.gz' \
  -minMaf 0.05 \
  -threads 16 \
  -iter 200 \
  -maf_iter 200 \
  -o $BASE_DIR'angsd/pcagnsd_bam_list_'$SAMPLE_SIZE'_uneven_coverage'
done
