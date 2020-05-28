#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/neutral_sim/}

BASE_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=12

## Get saf file
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -doSaf 1 \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL 2 \
      -P 2 \
      -doCounts 1 \
      -setMinDepth `awk "BEGIN {print $SAMPLE_SIZE*$COVERAGE}"` &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Get SFS from saf (this need to be broken up into two separate runs due to memory limiations)
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
      -P 2 \
      > $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Estimate theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -doThetas 1 \
      -doSaf 1 \
      -pest $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL 2 \
      -P 2 \
      -doCounts 1 \
      -setMinDepth `awk "BEGIN {print $SAMPLE_SIZE*$COVERAGE}"` &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Print per-SNP theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat print \
      $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
      > $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.tsv' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Do fixed window theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
      $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
      -win 10000 -step 10000 \
      -outnames $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Do per-chromosome average theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
      $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
      -outnames $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.idx' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
