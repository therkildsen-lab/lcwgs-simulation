#!/bin/bash
SAMPLE_SIZE=$1
COVERAGE=$2
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
## SNP calling
/workdir/programs/angsd0.931/angsd/angsd \
-b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
-anc $BASE_DIR'slim/ancestral.fasta' \
-out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 -doPost 1 -doVcf 1 \
-doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-P 4 -SNP_pval 1e-6 \
-setMinDepth 2 -minInd 1 -minQ 20 \
>& '/workdir/lcwgs-simulation/nohups/snp_calling_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.log'
## estimate SFS
/workdir/programs/angsd0.931/angsd/misc/realSFS $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
> $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs'
