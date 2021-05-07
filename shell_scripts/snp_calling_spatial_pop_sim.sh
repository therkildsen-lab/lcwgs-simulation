#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
GL_MODEL=${4:-1}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=6
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80}; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
    ## SNP calling
    /workdir/programs/angsd0.931/angsd/angsd \
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -out $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -GL $GL_MODEL -doGlf 2 -doMaf 1 -doMajorMinor 5 \
    -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -doIBS 2 -makematrix 1 -doCov 1 \
    -P 4 -SNP_pval 1e-6 -rmTriallelic 1e-6 \
    -setMinDepth 2 -minInd 1 -minMaf 0.05 -minQ 20 &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
