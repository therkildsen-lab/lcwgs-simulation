#!/bin/bash
BASE_DIR=${1:-/workdir/lcwgs-simulation/neutral_sim/rep_1/}
N_CORE_MAX=28

## SNP calling
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -out $BASE_DIR'angsd_gatk/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -GL 2 -doGlf 2 -doMaf 1 -doMajorMinor 5 \
    -doCounts 1 -doDepth 1 -dumpCounts 3 \
    -P 1 -SNP_pval 1e-6 -rmTriallelic 1e-6 \
    -setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
