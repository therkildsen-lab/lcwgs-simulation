#!/bin/bash
BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
REFERENCE=$2 # Path to reference genome
OUTDIR=$3 # Path to output files
OUTBASE=$4 # Basename for output files, e.g. GoM_Allsamples_MaxDepthMean1SD_MinQ20
SITES=$5
REGION=$6

/home/nt246_0001/Programs/angsd_0.912/angsd/angsd -b $BAMLIST -anc $REFERENCE -out $OUTDIR$OUTBASE -dosaf 1 -GL 1 -fold 1 -P 4 -sites $SITES -r $REGION >& $OUTDIR$OUTBASE'.log' 
