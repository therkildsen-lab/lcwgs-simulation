#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/}
SUBDIR=${3:-angsd}
GL_MODEL=${4:-1}
POP=${5:-1}

BASE_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=6

## Get saf file
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites' \
      -doSaf 1 \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL $GL_MODEL \
      -P 6 \
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
N_CORE_MAX=1
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.saf.idx' \
      -P 30 \
      -tole 1e-08 \
      -maxIter 1000 \
      > $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.sfs' &
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
N_CORE_MAX=6
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites' \
      -doThetas 1 \
      -doSaf 1 \
      -pest $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.sfs' \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL $GL_MODEL \
      -P 6 \
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
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \
      > $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.tsv' &
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
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \
      -win 10000 -step 10000 \
      -outnames $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.windowed_thetas.idx' &
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
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \
      -outnames $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.average_thetas.idx' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
