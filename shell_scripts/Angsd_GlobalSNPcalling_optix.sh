#!/bin/bash
BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
REFERENCE=$2 # Path to reference genome
OUTDIR=$3 # Path to output files
OUTBASE=$4 # Basename for output files, e.g. GoM_Allsamples_MaxDepthMean1SD_MinQ20
MINDP=$5 # 100
MAXDP=$6 # 1000
MININD=$7 # 100
MINQ=$8 #20
REGION=$9

/programs/angsd_20180926/angsd/angsd -b $BAMLIST -anc $REFERENCE -out $OUTDIR$OUTBASE -dosaf 1 -GL 1 -doGlf 3 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 16 -SNP_pval 1e-6 -setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD -minQ $MINQ -minMaf 0.05 -r $REGION  >& $OUTDIR$OUTBASE'.log'

#create a list of variant sites to use in downstream analyses 
cd $OUTDIR
gunzip -c $OUTBASE'.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > 'Global_SNPList_'$OUTBASE'.txt'
/home/nt246_0001/Programs/angsd_0.912/angsd/angsd sites index 'Global_SNPList_'$OUTBASE'.txt'

#also make it in regions format for downstream analyses
cut -f 1,2 'Global_SNPList_'$OUTBASE'.txt' | sed 's/\t/:/g' > 'Global_SNPList_'$OUTBASE'.regions'
