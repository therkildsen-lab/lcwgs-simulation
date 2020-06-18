#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=1
COUNT=0
# Generate SFS based on SAF			
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
		/workdir/programs/angsd0.931/angsd/misc/realSFS \
    $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    -P 32 -tole 3.16e-07 -maxIter 500 \
    > $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
# Estimate theta
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
      -bam $BASE_DIR'sample_lists/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -out $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -doThetas 1 \
      -doSaf 1 \
      -pest $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -GL 2 \
      -P 32 \
      -doCounts 1 \
      -setMinDepth `awk "BEGIN {print $SAMPLE_SIZE*$COVERAGE}"` &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
      wait
      COUNT=0
    fi
    done
done
wait
# Print per-SNP theta
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat print \
    $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    > $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.tsv'
  done
done
wait
# Do fixed window theta
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.windowed_thetas.idx'
  done
done
wait
## Do per-chromosome average theta
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -outnames $BASE_DIR'angsd_gatk/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.idx' 
  done
done
