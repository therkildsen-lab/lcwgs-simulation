#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
N_JOB_MAX=6
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    ## SNP calling
    /workdir/programs/angsd0.931/angsd/angsd \
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \
    -doCounts 1 -doDepth 1 -dumpCounts 1 \
    -P 6 -SNP_pval 1e-6 -rmTriallelic 1e-6 \
    -setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 \
    >& '/workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.log' &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
