#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/sim/rep_1/'
## SNP calling
/workdir/programs/angsd0.931/angsd/angsd \
-b $BASE_DIR'bam_list_50_1x.txt' \
-anc $BASE_DIR'ancestral_new.fasta' \
-out $BASE_DIR'bam_list_50_1x' \
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 -doPost 1 -doVcf 1 \
-doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-P 32 -SNP_pval 1e-6 \
-setMinDepth 10 -setMaxDepth 1000 -minInd 10 -minQ 20 -minMaf 0.01 \
>& /workdir/lcwgs-simulation/nohups/snp_calling_rep_1_bam_list_50_1x.log
## estimate SFS
/workdir/programs/angsd0.931/angsd/misc/realSFS $BASE_DIR'bam_list_50_1x.saf.idx' \
> $BASE_DIR'bam_list_50_1x.sfs'
