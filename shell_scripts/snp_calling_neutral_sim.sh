#!/bin/bash
SAMPLE_SIZE=$1
COVERAGE=$2
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
## SNP calling
/workdir/programs/angsd0.931/angsd/angsd \
-b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
-anc $BASE_DIR'slim/ancestral.fasta' \
-out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
-GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \
-doCounts 1 -doDepth 1 -dumpCounts 3 \
-P 1 -SNP_pval 1e-6 -rmTriallelic 1e-6 \
-setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 \
>& '/workdir/lcwgs-simulation/nohups/snp_calling_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.log'
