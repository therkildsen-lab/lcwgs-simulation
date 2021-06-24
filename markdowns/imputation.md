Evaluating the power of imputation without a reference panel to improve
genotype calls
================
Aryn Wilder

  - [Reformat bams for STITCH](#reformat-bams-for-stitch)
  - [Call SNPs and estimate MAFs in
    ANGSD](#call-snps-and-estimate-mafs-in-angsd)
  - [Define RunStitch.R](#define-runstitchr)
  - [Impute genotypes, genotype dosages and MAFS at SNPs in STITCH (may
    remove MAFS from this
    section)](#impute-genotypes-genotype-dosages-and-mafs-at-snps-in-stitch-may-remove-mafs-from-this-section)
  - [Impute genotypes and genotype dosages at SNPs in
    Beagle](#impute-genotypes-and-genotype-dosages-at-snps-in-beagle)
  - [Estimate ANGSD genotypes at called
    SNPs](#estimate-angsd-genotypes-at-called-snps)
  - [Get simulated (“true”) genotypes at SNPs from SLiM output
    fastas](#get-simulated-true-genotypes-at-snps-from-slim-output-fastas)
  - [Convert simulated genotypes to
    dosages](#convert-simulated-genotypes-to-dosages)
  - [Get dosage from posterior probabilities in
    ANGSD](#get-dosage-from-posterior-probabilities-in-angsd)
  - [Convert ANGSD posterior GPs to
    dosage](#convert-angsd-posterior-gps-to-dosage)
  - [Estimate MAFS from Beagle GPs](#estimate-mafs-from-beagle-gps)
  - [Estimate MAFS from STITCH GPs](#estimate-mafs-from-stitch-gps)
  - [Get MAFs from SLiM](#get-mafs-from-slim)
  - [Compare genotype dosage estimated in ANGSD, STITCH and Beagle
    compared to true dosage simulated in
    SLiM](#compare-genotype-dosage-estimated-in-angsd-stitch-and-beagle-compared-to-true-dosage-simulated-in-slim)
  - [Compare genotypes estimated in ANGSD, STITCH and Beagle compared to
    true GTs simulated in
    SLiM](#compare-genotypes-estimated-in-angsd-stitch-and-beagle-compared-to-true-gts-simulated-in-slim)
  - [Compile and plot results from all three simulated
    recombination/mutation
    scenarios](#compile-and-plot-results-from-all-three-simulated-recombinationmutation-scenarios)
      - [Compare MAFs estimated in ANGSD, STITCH and Beagle to those
        simulated in
        SLiM](#compare-mafs-estimated-in-angsd-stitch-and-beagle-to-those-simulated-in-slim)
      - [Plot r2 of genotypes estimated in ANGSD, STITCH and Beagle
        compared to true genotypes simulated in
        SLiM](#plot-r2-of-genotypes-estimated-in-angsd-stitch-and-beagle-compared-to-true-genotypes-simulated-in-slim)

All of the steps for testing imputation were done on a different server
than the one used for the simulations. The following files were
transferred:

1.  bam files in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_high_r/rep_1/bam/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_hi/bams/`  
    `ancestral.fasta` in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_high_r/rep_1/` were
    moved to `~/USS/aryn/LCWGS/LCWGS_hi/`  
    fasta files in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_high_r/rep_1/fasta/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_hi/fastas/`

2.  bam files in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_medium_r/rep_1/bam/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_med/bams/`  
    `ancestral.fasta` in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_medium_r/rep_1/` were
    moved to `~/USS/aryn/LCWGS/LCWGS_med/`  
    fasta files in
    `/workdir/lcwgs-simulation/neutral_sim_high_mu_medium_r/rep_1/fasta/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_med/fastas/`

3.  bam files in
    `/workdir/lcwgs-simulation/neutral_sim_medium_mu_low_r/rep_1/bam/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_low/bams/`  
    ancestral.fasta in
    `/workdir/lcwgs-simulation/neutral_sim_medium_mu_low_r/rep_1/` were
    moved to `~/USS/aryn/LCWGS/LCWGS_low/`  
    fasta files in
    `/workdir/lcwgs-simulation/neutral_sim_medium_mu_low_r/rep_1/fasta/`
    were moved to `~/USS/aryn/LCWGS/LCWGS_low/fastas/`

First run each step of the pipeline for the three simulated
recombination/mutation scenarions (LCWGS\_low, LCWGS\_med, and
LCWGS\_hi)

## Reformat bams for STITCH

``` bash
#!/bin/bash

SCENARIO=$1

cd '~/USS/aryn/LCWGS/'$SCENARIO'/bams'

for COV in 1 2 4; do
  for BAMFILE in `ls *'_sorted_'$COV'x.bam'`; do
    FILENAME=`echo ${BAMFILE##*/} | sed 's/.bam//g'`
    ./reformat.sh in=$BAMFILE out=$FILENAME"_rf.bam" sam=1.3
    samtools index $FILENAME"_rf.bam"
  done

  for N in 25 100 250 500 1000; do
    find $(pwd) -maxdepth 1 -type f -name '*_sorted_'$COV'x_rf.bam' | shuf | head -n $N > \
    '~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
  done
done

cd '~/USS/aryn/LCWGS/'$SCENARIO
```

Identify variant sites for imputation

## Call SNPs and estimate MAFs in ANGSD

``` bash
#!/bin/bash
SCENARIO=$1

samtools faidx ancestral.fasta 
mkdir '~/USS/aryn/LCWGS/'$SCENARIO'/angsd'

for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    BAMLIST='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
    POSFILE='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_pos.txt'
    
    #Call SNPs in angsd
    angsd -b $BAMLIST \
    -anc '~/USS/aryn/LCWGS/'$SCENARIO'/ancestral.fasta' \
    -out '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x' -GL 1 \
    -doGlf 2 -doMaf 1 -doMajorMinor 5 -doCounts 1 -doDepth 1 -dumpCounts 3 -P 12 \
    -SNP_pval 1e-6 -rmTriallelic 1e-6 -setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20
    
    #create list of SNP positions
    zcat '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x.mafs.gz' \
    | tail -n +2 | awk '{print $1"\t"$2"\t"$3"\t"$4}' > $POSFILE
  done
done
```

## Define RunStitch.R

``` r
shell_script <- "args=commandArgs(trailingOnly = TRUE)
N=args[1]
COV=args[2]
K=as.numeric(args[3])
NGEN=10
S=4
BAMLIST=args[4]
POSFILE=args[5]
SAMPLENAMES=args[6]
OUTDIR=args[7]
dir=getwd()


library('STITCH')

STITCH(chr='rep_1', bamlist=BAMLIST, posfile=POSFILE, outputdir=OUTDIR, K=K, nGen=NGEN, nCores=14, S=S, sampleNames_file=SAMPLENAMES, 
tempdir = 'tmpfiles')"
write_lines(shell_script, "RunStitch.R")
```

## Impute genotypes, genotype dosages and MAFS at SNPs in STITCH (may remove MAFS from this section)

``` bash

mkdir '~/USS/aryn/LCWGS/'$SCENARIO'/stitch'
mkdir '~/USS/aryn/LCWGS/'$SCENARIO'/tmpfiles'

SCENARIO=$1
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    K=30
    NGEN=10
    S=4
    BAMLIST='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
    POSFILE='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_pos.txt'
    SAMPLENAMES='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_samplelist.txt'
    OUTDIR='~/USS/aryn/LCWGS/'$SCENARIO'/stitch/stitch_k'$K'_ngen'$NGEN'_s'$S'_n'$N'_'$COV'x'

    #Create samplelist
    if [ ! -f "$SAMPLENAMES" ]
    then
      for BAM in `cat $BAMLIST`; do
        SAMPLE=`echo ${BAM##*/} | cut -d_ -f1-2`
        echo $SAMPLE >> $SAMPLENAMES
      done
    fi
    
    Rscript RunStitch.R $N $COV $K $BAMLIST $POSFILE $SAMPLENAMES $OUTDIR
    
    #get genotypes from STITCH for comparing with simulated genotypes
    tabix $OUTDIR'/stitch.rep_1.vcf.gz'
    zcat $OUTDIR'/stitch.rep_1.vcf.gz' | head -n 15 | grep "#CHROM" | cut -f1,2,10- > $OUTDIR'/stitch.rep_1.GT'
    bcftools query -f '%CHROM\t%POS\t[\t%TGT]\n' $OUTDIR'/stitch.rep_1.vcf.gz' >> $OUTDIR'/stitch.rep_1.GT'
    sed -e 's.C\/A.AC.g' -e 's.G\/A.AG.g' -e 's.T\/A.AT.g' -e 's.G\/C.CG.g' -e 's.T\/C.CT.g' -e 's.T\/G.GT.g' -e 's.\/..g' $OUTDIR'/stitch.rep_1.GT' > $OUTDIR'/stitch_k'$K'_ngen'$NGEN'_s'$S'_n'$N'_'$COV'x_ordered.geno'
    gzip $OUTDIR'/stitch.rep_1.GT' &
    
    #get genotype dosage from STITCH for comparing with simulated genotypes
    vcftools --gzvcf $OUTDIR'/stitch.rep_1.vcf.gz' \
    --extract-FORMAT-info DS \
    --out $OUTDIR'/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x.dosage'

  done
done
```

## Impute genotypes and genotype dosages at SNPs in Beagle

``` bash
SCENARIO=$1

mkdir '~/USS/aryn/LCWGS/'$SCENARIO'/beagle'

for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    #run imputation
    java -Djava.io.tmpdir='tmpfiles' -jar beagle.jar \
    like='~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gz' \
    out='~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed'
    
    #call genotypes with posterior prob>0.9
    zcat '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gz.gprobs.gz' | \
    java -jar gprobs2beagle.jar 0.90 N > '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gt'
    
    #change numbered alleles to bases for comparing with simulated genotypes
    tail -n +2 '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gt' | sed -e 's.0.A.g' -e 's.1.C.g' -e 's.2.G.g' -e 's.3.T.g' \
    | awk '{for (i=3;i<=NF;i+=2) {printf ($i""$(i+1)"\t")} {printf "\n"}}' \
    > '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.geno'
    
    #alphabetize alleles in genotypes for comparing with simulated genotypes
    sed -e 's.CA.AC.g' -e 's.GA.AG.g' -e 's.TA.AT.g' -e 's.GC.CG.g' -e 's.TC.CT.g' -e 's.TG.GT.g' \
    '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.geno' > \
    '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle_ordered.geno'
    paste <(tail -n +2 '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gt' | cut -d" " -f1-2) \
    '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle_ordered.geno' > \
    '~/USS/aryn/LCWGS/'$SCENARIO'/beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle_ordered.geno2'
  done
done
```

## Estimate ANGSD genotypes at called SNPs

``` bash
SCENARIO=$1
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    BAMLIST='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
    POSFILE='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_pos.txt'
    
    angsd sites index $POSFILE
    
    #call genotypes with posterior prob>0.9 at all SNPs
    angsd -b $BAMLIST \
    -anc '~/USS/aryn/LCWGS/'$SCENARIO'/ancestral.fasta' -postCutoff 0.90 -doPost 1 -doMaf 1 -GL 1 -dogeno 5 \
    -out '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior' \
    -doMajorMinor 3 -sites $POSFILE -P 12 

    #reformat to compare with "true" genotypes
    zcat '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno.gz' | \
    cut -f1-2 > '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/markernames'
    zcat '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno.gz' | \
    cut -f 3- > '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno'
    paste '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/markernames' \
    '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno' | gzip \
    > '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/Angsd_neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno.gz'
    rm '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno' \
    '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/markernames'

    #order alleles of genotypes alphabetically
    zcat '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/Angsd_neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior.geno.gz' | \
    sed -e 's.CA.AC.g' -e 's.GA.AG.g' -e 's.TA.AT.g' -e 's.GC.CG.g' -e 's.TC.CT.g' -e 's.TG.GT.g' \
    > '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/Angsd_neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior_ordered.geno'
  done
done
```

## Get simulated (“true”) genotypes at SNPs from SLiM output fastas

``` bash
SCENARIO=$1
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    BAMLIST='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
    POSFILE='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_pos.txt'

    REGFILE=`echo $POSFILE | sed 's/txt/regions/g'`
    awk '{print $1":"$2"-"$2}' $POSFILE > $REGFILE

    mkdir '~/USS/aryn/LCWGS/'$SCENARIO'/True_N'$N'_'$COV'x'
    OUTDIR='~/USS/aryn/LCWGS/'$SCENARIO'/True_N'$N'_'$COV'x/'
    
    awk '{print $1"\t"$2}' $POSFILE > $OUTDIR'True_N'$N'.geno'
    for BAM in `cat $BAMLIST`; do
      SAMPLE=`echo ${BAM##*/} | cut -d_ -f1-2`
      echo $SAMPLE
      FASTAS=`samtools view -H $BAM | grep "\-sam" | cut -f 4 | cut -d" " -f6`
      FASTA1=`echo $FASTAS | cut -d" " -f1`
      FASTA1=`echo ${FASTA1##*/}`
      FASTA2=`echo $FASTAS | cut -d" " -f2`
      FASTA2=`echo ${FASTA2##*/}`
    
      samtools faidx '~/USS/aryn/LCWGS/'$SCENARIO'/fastas/'$FASTA1'.gz' -r $REGFILE | grep -v ">" > \
      $OUTDIR$SAMPLE".gt1"
      samtools faidx '~/USS/aryn/LCWGS/'$SCENARIO'/fastas/'$FASTA2'.gz' -r $REGFILE | grep -v ">" > \
      $OUTDIR$SAMPLE".gt2"
    
      paste -d'\0' $OUTDIR$SAMPLE".gt1" \
      $OUTDIR$SAMPLE".gt2" > \
      $OUTDIR$SAMPLE".gt"
      paste $OUTDIR'True_N'$N'.geno' \
      $OUTDIR$SAMPLE'.gt' > \
      $OUTDIR'True_N'$N'.geno1'
      mv $OUTDIR'True_N'$N'.geno1' \
      $OUTDIR'True_N'$N'.geno'
      rm $OUTDIR$SAMPLE".gt"*
    done

    sed -e 's.CA.AC.g' -e 's.GA.AG.g' -e 's.TA.AT.g' -e 's.GC.CG.g' -e 's.TC.CT.g' -e 's.TG.GT.g' \
    $OUTDIR'True_N'$N'.geno' > \
    $OUTDIR'True_N'$N'_'$COV'x_ordered.geno'
  done
done
```

## Convert simulated genotypes to dosages

``` r
library(stringr)
args=commandArgs(trailingOnly = TRUE)
scenario=args[1]
setwd(paste0('~/USS/aryn/LCWGS/',scenario,collapse=''))

for (x in c(1,2,4)){
  for (N in c(25,100,250,500,1000)){

    slim=read.table(paste(c('True_N',N,'_',x,'x/True_N',N,'_',x,'x_ordered.geno'),collapse=''),colClasses='character')
    mafs=read.table(paste(c('angsd/neutral_sim_with_replacement_n',N,'_',x,'x.mafs.gz'),collapse=''),header=T)
    row.names(slim)=slim$V2
    slim=slim[as.character(mafs$position),]
    
    slimgt=matrix(nrow=nrow(mafs),ncol=ncol(slim))
    slimgt[,1]=mafs[,1]
    slimgt[,2]=mafs[,2]
    
    for (i in 1:nrow(mafs)){
        anc=mafs$anc[i]
        slimgt[i,3:ncol(slimgt)] =sapply(slim[i,3:ncol(slim)],function(x) 2-(str_count(as.character(x), as.character(anc))))
    }
    slimgt=as.data.frame(slimgt)
    
    write.table(slimgt,paste(c('True_N',N,'_',x,'x/True_N',N,'_',x,'x.dosage'),collapse=''),sep='\t',row.names=F,col.names=F,quote=F)
  }
}
```

## Get dosage from posterior probabilities in ANGSD

first estimate posterior GPs

``` bash
SCENARIO=$1
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    BAMLIST='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_rf_bamlist.txt'
    POSFILE='~/USS/aryn/LCWGS/'$SCENARIO'/neutral_sim_with_replacement_n'$N'_'$COV'x_pos.txt'

    angsd -b $BAMLIST \
    -anc '~/USS/aryn/LCWGS/'$SCENARIO'/ancestral.fasta' -doPost 1 -doMaf 1 -GL 1 -dogeno 8 \
    -out '~/USS/aryn/LCWGS/'$SCENARIO'/angsd/neutral_sim_with_replacement_n'$N'_'$COV'x_AFprior_dopost8' \
    -doMajorMinor 3 -sites $POSFILE -P 12 >& 'angsd_dopost_n'$N'_'$COV'x.log'
  done
done
```

## Convert ANGSD posterior GPs to dosage

``` r
args=commandArgs(trailingOnly = TRUE)
scenario=args[1]
setwd(paste0('~/USS/aryn/LCWGS/',scenario,collapse=''))

for (cov in c(1,2,4)){
  for (n in c(25,100,250,500,1000)){
    chr=read.table(paste(c('angsd/neutral_sim_with_replacement_n',n,'_',cov,'x_AFprior_dopost8.geno.gz'),collapse=''))
    
    chrgt=as.data.frame(matrix(ncol=(ncol(chr)-2)/3+2,nrow=nrow(chr)))
    chrgt[,1:2]=chr[,1:2]
    for (i in 1:nrow(chr)){
        for (j in 1:((ncol(chr)-2)/3)){
                chrgt[i,j+2]=chr[i,j*3]*0+chr[i,1+j*3]*1+chr[i,2+j*3]*2
        }
    }

    write.table(chrgt,paste(c('angsd/neutral_sim_with_replacement_n',n,'_',cov,'x_AFprior_dosage.txt'),collapse=''),col.names=F,row.names=F,sep='\t',quote=F)
  }
}
```

## Estimate MAFS from Beagle GPs

``` bash
SCENARIO=$1
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do

    angsd -beagle 'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gz.gprobs.gz' \
    -doMaf 4 -P 6 -intName 0 -fai ancestral.fasta.fai \
    -out 'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x_beagle' \
    &> 'beaglemaf_'$N'_'$COV'.log'
    
    zcat 'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x.beagle.gz.gprobs.gz' | cut -d" " -f1 | cut -d_ -f3 > \
    'beagle/positions_'$N'_'$COV
    paste 'beagle/positions_'$N'_'$COV <(zcat 'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x_beagle.mafs.gz' | cut -f3-) > \
    'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x_beagle.mafs.txt'
    rm 'beagle/positions_'$N'_'$COV 'beagle/imputed.neutral_sim_with_replacement_n'$N'_'$COV'x_beagle.mafs.gz'
  done
done
```

## Estimate MAFS from STITCH GPs

``` bash
for COV in 1 2 4; do
  for N in 25 100 250 500 1000; do
    K=30
    BASENAME='stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x'

    vcftools --gzvcf 'stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/stitch.rep_1.vcf.gz' \
    --extract-FORMAT-info GP \
    --out $BASENAME'.gprob'

    cut -f 3- $BASENAME'.gprob.GP.FORMAT' | sed 's.\,.\t.g' > $BASENAME'.gprob.GPs'
    awk '{print $1"_"$2}' $BASENAME'.gprob.GP.FORMAT' | tail -n +2 > $BASENAME'.gprob.GP.sites'
    zcat 'stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/stitch.rep_1.vcf.gz' | grep -v "#" | cut -f4-5 \
    > 'stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/alleles'
    paste $BASENAME'.gprob.GP.sites' 'stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/alleles' \
    > $BASENAME'.gprob.GP.sitesalleles'
    paste $BASENAME'.gprob.GP.sitesalleles' <(tail -n +2 $BASENAME'.gprob.GPs') | gzip > $BASENAME'.GP.gz'
    zcat $BASENAME'.GP.gz' | head -n 1 > $BASENAME'temp'
    zcat $BASENAME'.GP.gz' >> $BASENAME'temp'
    gzip $BASENAME'temp'
    mv $BASENAME'temp.gz' $BASENAME'.GP.gz'
    rm $BASENAME'.gprob.GPs' $BASENAME'.gprob.GP.sites' $BASENAME'.gprob.GP.sitesalleles' \
    $BASENAME'.gprob.GP.FORMAT' 'stitch/stitch_k'$K'_ngen10_s4_n'$N'_'$COV'x/alleles'
    
    angsd -doMaf 4 -P 6 -intName 0 \
    -beagle $BASENAME'.GP.gz' \
    -fai ancestral.fasta.fai \
    -out $BASENAME'.GP'

    zcat $BASENAME'.GP.gz' | cut -f1 | cut -d_ -f3 | tail -n +2 > $BASENAME'positions'
    paste $BASENAME'positions' <(zcat $BASENAME'.GP.mafs.gz' | tail -n +2 | cut -f3-) > $BASENAME'.mafs.txt'
    rm $BASENAME'positions' 
  done
done
```

## Get MAFs from SLiM

``` r
library(tidyverse)
args=commandArgs(trailingOnly = TRUE)
scenario=args[1]
setwd(paste0('~/USS/aryn/LCWGS/',scenario,collapse=''))

ancestral <- read_csv("ancestral.fasta")[[1]] %>%
  str_split(pattern="") %>%
  .[[1]] %>%
  bind_cols(ancestral=., position=1:30000000)

## Read in the mutation file outputted by SLiM
mutations <- read_delim("mutations.txt", delim = " ", col_names = F) %>%
  transmute(type=X6, position=X7+1, base=X13, frequency=X12/2000) %>%
  left_join(ancestral, by="position") %>%
  group_by(type, position, ancestral, base) %>%
  summarise(frequency=sum(frequency)) %>%
  ungroup()
## Read in the substitutions file outputted by SLiM
substitutions <- read_delim("substitutions.txt", delim = " ", skip=2, col_names = F) %>%
  transmute(type=X3, position=X4+1, base=X10, generation=X9) %>%
  group_by(type, position) %>%
  filter(generation==max(generation)) %>%
  ungroup() %>%
  left_join(ancestral, by="position") %>%
  select(-generation) %>%
  filter(base!=ancestral) %>%
  arrange(position)
## Join mutations and substitutions in a temp table
mutations_final_temp <-  mutations %>%
  spread(key = base, value=frequency) %>%
  full_join(substitutions, by=c("position", "type", "ancestral")) %>%
  arrange(position) %>%
  mutate(base=ifelse(is.na(base), ancestral, base)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(frequency=1-`A` -`C` -`G` -`T`)
## More wrangling
mutations_final <- mutations_final_temp[1:7] %>%
  gather(key=base, value=frequency, 4:7) %>%
  bind_rows(mutations_final_temp[c(1:3, 8:9)]) %>%
  mutate(frequency=ifelse(base==ancestral, 0, frequency)) %>%
  group_by(type, position, ancestral) %>%
  filter(frequency!=0) %>%
  summarise(frequency=sum(frequency), base=paste0(base, collapse = "")) %>%
  ungroup() %>%
  filter(frequency!=1)

write.table(mutations_final,'slim_mafs.txt',quote=F,row.names=F,sep='\t')
```

## Compare genotype dosage estimated in ANGSD, STITCH and Beagle compared to true dosage simulated in SLiM

``` r
library(stringr)
args=commandArgs(trailingOnly = TRUE)
scenario=args[1]
setwd(paste0('~/USS/aryn/LCWGS/',scenario,collapse=''))

for (inds in c(25,100,250,500,1000)){
    for (cov in c(1,2,4)){
    K=30
    
    stitch=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/stitch/stitch_k',K,'_ngen10_s4_n',inds,'_',cov,'x/stitch_k',K,'_ngen10_s4_n',inds,'_',cov,'x.dosage.DS.FORMAT'),collapse=''),colClasses='character',na.strings='..',header=T)
    slim=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/True_N',inds,'_',cov,'x/True_N',inds,'_',cov,'x.dosage'),collapse=''))#,colClasses='character'
    mafs=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/angsd/neutral_sim_with_replacement_n',inds,'_',cov,'x.mafs.gz'),collapse=''),header=T)
    bgl=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/beagle/imputed.neutral_sim_with_replacement_n',inds,'_',cov,'x.beagle.gz.dose.gz'),collapse=''),header=T,nrow=300000)
    row.names(bgl)=str_split_fixed(bgl$marker,"_",3)[,3]
    ang=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/angsd/neutral_sim_with_replacement_n',inds,'_',cov,'x_AFprior_dosage.txt'),collapse=''),header=F)
    mafs=merge(mafs,stitch[,1:2],by.x=c('chromo','position'),by.y=c('CHROM','POS'))
    mafs=merge(mafs,slim[,1:2],by.x='position',by.y='V2')
    mafs=mafs[order(mafs$position),]
    mafs$anc=as.character(mafs$anc)
    row.names(slim)=slim$V2
    row.names(stitch)= stitch$POS
    row.names(ang)=ang$V2
    slim=slim[as.character(mafs$position),]
    stitch=stitch[as.character(mafs$position),]
    bgl=bgl[as.character(mafs$position),]
    ang=ang[as.character(mafs$position),]
    
        
    impr2=matrix(nrow=20,ncol=5)
    impr2[,1]=seq(0,0.95,0.05)
    r=0

    #r2 between stitch and simulated
for (j in seq(0,0.95,0.05)){
    r=r+1
    idx=which(mafs$knownEM>=j & mafs$knownEM<(j+0.05))
    r2=matrix(nrow=length(idx),ncol=1)
    for (k in 1:length(idx)){
    r2[k]=summary(lm(as.numeric(unlist(stitch[idx[k],3:ncol(stitch)]))~as.numeric(unlist(slim[idx[k],3:ncol(slim)]))))$r.squared
    }
impr2[r,2]=sum(!is.na(r2))
impr2[r,3]=round(mean(r2,na.rm=T),6)
}
    
    r=0
    #r2 between beagle and simulated
    for (j in seq(0,0.95,0.05)){
        r=r+1
        idx=which(mafs$knownEM>=j & mafs$knownEM<(j+0.05))
        r2=matrix(nrow=length(idx),ncol=1)
        for (k in 1:length(idx)){
        r2[k]=summary(lm(unlist(bgl[idx[k],4:ncol(bgl)])~unlist(slim[idx[k],3:ncol(slim)])))$r.squared
        }
    impr2[r,4]=round(mean(r2,na.rm=T),6)
    }
    
    r=0
    #r2 between angsd and simulated
    for (j in seq(0,0.95,0.05)){
            r=r+1
        idx=which(mafs$knownEM>=j & mafs$knownEM<(j+0.05))
            r2=matrix(nrow=length(idx),ncol=1)
            for (k in 1:length(idx)){
            r2[k]=summary(lm(unlist(ang[idx[k],3:ncol(ang)])~unlist(slim[idx[k],3:ncol(slim)])))$r.squared
            }
    impr2[r,5]=round(mean(r2,na.rm=T),6)
    }
    
    impr2=as.data.frame(impr2)
    names(impr2)=c('maf','sites','stitch_r2','beagle_r2','angsd_r2')
    write.table(impr2,paste(c('~/USS/aryn/LCWGS/',scenario,'/DosageR2_N',inds,'_',cov,'x.txt'),collapse=''),sep='\t',row.names=F,quote=F)
  }
}
```

## Compare genotypes estimated in ANGSD, STITCH and Beagle compared to true GTs simulated in SLiM

``` r
args=commandArgs(trailingOnly = TRUE)
scenario=args[1]
setwd(paste0('~/USS/aryn/LCWGS/',scenario,collapse=''))

mutations_final=read.table('slim_mafs.txt',header=T,sep='\t')
hets=c('AC','AG','AT','CG','CT','GT')
homs=c('AA','CC','GG','TT')

evalGT=matrix(nrow=12,ncol=8)
i=0
for (inds in c(25,100,250,500,1000)){
    for (cov in c(1,2,4)){
    i=i+1
    k=30
    
    simgt=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/True_N',inds,'_',cov,'x/True_N',inds,'_',cov,'x_ordered.geno'),collapse=''),sep='\t',colClasses='character')
    row.names(simgt)=simgt[,2]
    impgt=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/beagle/imputed.neutral_sim_with_replacement_n',inds,'_',cov,'x.beagle_ordered.geno'),collapse=''),sep='\t',colClasses='character',na.strings='NN')
    # impgt=apply(impgt,2,as.character)
    impgt=impgt[,c(1:(ncol(impgt)-1))]#2,4
    row.names(impgt)=simgt[,2]
    anggt=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/angsd/Angsd_neutral_sim_with_replacement_n',inds,'_',cov,'x_AFprior_ordered.geno'),collapse=''),sep='\t',colClasses=rep('character',inds
    +2),na.strings='NN')
    anggt=anggt[,c(1,2,5:(ncol(anggt)-1))]
    row.names(anggt)= anggt[,2]
    stgt=read.table(paste(c('~/USS/aryn/LCWGS/',scenario,'/stitch/stitch_k',k,'_ngen10_s4_n',inds,'_',cov,'x/stitch_k',k,'_ngen10_s4_n',inds,'_',cov,'x_ordered.geno'),collapse=''),sep='\t',colClasses='character',na.string='..')
    stgt= stgt[,c(1:2,4:ncol(stgt))]
    row.names(stgt)= stgt[,2]
    
    missing=simgt$V2[which(simgt$V2 %in% anggt$V2==F)]
    missmat=matrix(nrow=length(missing),ncol=ncol(anggt))
    missmat[,2]=missing
    missmat[,1]='rep_1'
    missmat=as.data.frame(missmat)
    names(missmat)=names(anggt)
    anggt=rbind(anggt,missmat)
    anggt= anggt[order(as.numeric(anggt$V2)),]
    
    angcorr=sum(simgt[,3:ncol(simgt)]==anggt[,3:ncol(anggt)],na.rm=T)/sum(!is.na(anggt[,3:ncol(anggt)]))
    impcorr=sum(simgt[,3:ncol(simgt)]==impgt,na.rm=T)/sum(!is.na(impgt))#[,3:ncol(impgt)]
    stcorr=sum(simgt[,3:ncol(simgt)]==stgt[,3:ncol(stgt)],na.rm=T)/sum(!is.na(stgt[,3:ncol(stgt)]))
    angcalled=sum(!is.na(anggt[,3:ncol(anggt)]))/(sum(is.na(anggt[,3:ncol(anggt)]))+sum(!is.na(anggt[,3:ncol(anggt)])))
    impcalled=sum(!is.na(impgt))/(sum(is.na(impgt))+sum(!is.na(impgt)))#[,3:ncol(impgt)]
    stcalled=sum(!is.na(stgt))/(sum(is.na(stgt))+sum(!is.na(stgt)))#[,3:ncol(impgt)]
    
    evalGT[i,1:2]=c(inds,cov)
    evalGT[i,3:8]=c(angcorr,impcorr,stcorr,angcalled,impcalled,stcalled)
    
    #evaluate by MAF bins
    Corrcalled=matrix(ncol=11,nrow=20)
    HetHom=matrix(nrow=20,ncol=8)
    for (j in seq(0,0.95,0.05)){
    hifrq=mutations_final$position[which(mutations_final$frequency>=j & mutations_final$frequency<(j+0.05))] #see ValidateImputations.R to get this file
    hifrq=as.character(hifrq)
    hifrq=hifrq[hifrq %in% row.names(simgt)] 
    
    Corrcalled[(j+0.05)*20,1]=sum(simgt[hifrq,3:ncol(simgt)]==impgt[hifrq,],na.rm=T)/sum(!is.na(impgt[hifrq,]))
    Corrcalled[(j+0.05)*20,2]=sum(!is.na(impgt[hifrq,]))/(sum(!is.na(impgt[hifrq,]))+sum(is.na(impgt[hifrq,])))
    Corrcalled[(j+0.05)*20,3]=sum(simgt[hifrq[hifrq %in% row.names(anggt)],3:ncol(simgt)]==anggt[hifrq[hifrq %in% row.names(anggt)],3:ncol(anggt)],na.rm=T)/sum(!is.na(anggt[hifrq[hifrq %in% row.names(anggt)],3:ncol(anggt)]))
    Corrcalled[(j+0.05)*20,4]=sum(!is.na(anggt[hifrq[hifrq %in% row.names(anggt)],3:ncol(anggt)]))/(sum(!is.na(anggt[hifrq[hifrq %in% row.names(anggt)],3:ncol(anggt)]))+sum(is.na(anggt[hifrq[hifrq %in% row.names(anggt)],3:ncol(anggt)])))
    Corrcalled[(j+0.05)*20,5]=sum(simgt[hifrq,3:ncol(simgt)]==stgt[hifrq,3:ncol(stgt)],na.rm=T)/sum(!is.na(stgt[hifrq,3:ncol(stgt)]))
    Corrcalled[(j+0.05)*20,6]=sum(!is.na(stgt[hifrq,3:ncol(stgt)]))/(sum(!is.na(stgt[hifrq,3:ncol(stgt)]))+sum(is.na(stgt[hifrq,3:ncol(stgt)])))
    Corrcalled[(j+0.05)*20,7]=length(hifrq)
    
    hifrq=mutations_final$position[which(mutations_final$frequency>=j & mutations_final$frequency<(j+0.05))] #see ValidateImputations.R to get this file
    hifrq=as.character(hifrq)
    hifrq=hifrq[hifrq %in% simgt[,2]] 
    
    Corrcalled[(j+0.05)*20,8]=sum(rowSums(apply(simgt[hifrq,3:ncol(simgt)],2, function(x) x %in% hets),na.rm=T))/sum(!is.na(simgt[hifrq,3:ncol(simgt)]))
    Corrcalled[(j+0.05)*20,9]=sum(rowSums(apply(anggt[hifrq,3:ncol(anggt)],2, function(x) x %in% hets),na.rm=T))/sum(!is.na(anggt[hifrq,3:ncol(anggt)]))
    Corrcalled[(j+0.05)*20,10]=sum(rowSums(apply(impgt[hifrq,],2, function(x) x %in% hets),na.rm=T))/sum(!is.na(impgt[hifrq,]))
    Corrcalled[(j+0.05)*20,11]=sum(rowSums(apply(stgt[hifrq,3:ncol(stgt)],2, function(x) x %in% hets),na.rm=T))/sum(!is.na(stgt[hifrq,3:ncol(stgt)]))
    
    }
    Corrcalled=as.data.frame(Corrcalled)
    names(Corrcalled)=c('Beagle_corr','Beagle_called','Angsd_corr','Angsd_called','Stitch_corr','Stitch_called','sites','True_het','Angsd_het','Beagle_het','Stitch_het')
    write.table(Corrcalled,paste(c('CalledGT_byMAF_N',inds,'_',cov,'x.txt'),collapse=''),quote=F,row.names=F,sep='\t')
  
  }
}
evalGT =as.data.frame(evalGT)
names(evalGT)=c('N','x','Angsd_corr','Beagle_corr','Stitch_corr','Angsd_called','Beagle_called','Stitch_called')
```

## Compile and plot results from all three simulated recombination/mutation scenarios

### Compare MAFs estimated in ANGSD, STITCH and Beagle to those simulated in SLiM

``` r
library(Metrics)
setwd('~/USS/aryn/LCWGS')

i=0
MAFimputation=matrix(ncol=9,nrow=45)
for (scenario in c('LCWGS_low','LCWGS_med','LCWGS_hi')){
  for (cov in c(1,2,4)){
    for (n in c(25,100,250,500,1000)){
i=i+1
k=30

MAFimputation[i,1]=scenario
MAFimputation[i,2]=cov
MAFimputation[i,3]=n

sti=read.table(paste(c(scenario,'/stitch/stitch_k',k,'_ngen10_s4_n',n,'_',cov,'x/stitch_k',k,'_ngen10_s4_n',n,'_',cov,'x.mafs.txt'),collapse=''))
names(sti)=c('position','stitch_major','stitch_minor','stitch_maf','stitch_nInd')
ang=read.table(paste(c(scenario,'/angsd/neutral_sim_with_replacement_n',n,'_',cov,'x.mafs.gz'),collapse=''),header=T)
names(ang)[2:8]=c('position','ang_major','ang_minor','ang_ancestral','ang_maf','ang_pK.EM','ang_nInd')
tru=read.table(paste(c(scenario,'/slim_mafs.txt'),collapse=''),header=T)
bgl=read.table(paste(c(scenario,'/beagle/imputed.neutral_sim_with_replacement_n',n,'_',cov,'x_beagle.mafs.txt'),collapse=''),header=T)
names(bgl)=c('position','bgl_major','bgl_minor','bgl_maf','bgl_nInd')
tru=merge(tru,sti,by='position')
tru=merge(tru,ang[,c(2:6,8)],by='position')
tru=merge(tru,bgl,by='position')

MAFimputation[i,4]=summary(lm(tru$frequency~tru$ang_maf))$r.squared
MAFimputation[i,5]=summary(lm(tru$frequency~tru$stitch_maf))$r.squared
MAFimputation[i,6]=summary(lm(tru$frequency~tru$bgl_maf))$r.squared

MAFimputation[i,7]=rmse(tru$frequency,tru$ang_maf)
MAFimputation[i,8]=rmse(tru$frequency,tru$stitch_maf)
MAFimputation[i,9]=rmse(tru$frequency,tru$bgl_maf)
    }}}
MAFimputation=as.data.frame(MAFimputation)
names(MAFimputation)=c('Scenario','Cov','N','Angsd_r2','Stitch_r2','Beagle_r2','Angsd_RMSE','Stitch_RMSE','Beagle_RMSE')
MAFimputation[4:9]=apply(MAFimputation[4:9],2,function(x) as.numeric(as.character(x)))
write.table(MAFimputation,'MAFimputationImprovement.txt',sep='\t',row.names=F,quote=F)
```

### Plot r2 of genotypes estimated in ANGSD, STITCH and Beagle compared to true genotypes simulated in SLiM

``` r
library(ggplot2)
library(reshape)
library(stringr)
library(lemon)

#################################### r2 by MAF ###################################
DS=NULL
for (scenario in c('LCWGS_low','LCWGS_med','LCWGS_hi')){
  setwd(paste(c('~/USS/aryn/LCWGS/',scenario),collapse=''))
    for (N in c(25,100,250,500,1000)){
      for (x in c(1,2,4)){
      sDS=read.table(paste(c('DosageR2_N',N,'_',x,'x.txt'),collapse=''),header=T)
      mDS=matrix(nrow=10,ncol=7)
      for (i in 1:10){
        k=20-(i-1)
        for (j in 3:5){
            mDS[i,j]=weighted.mean(c(sDS[i,j],sDS[k,j]),c(sDS[i,2],sDS[k,2]))
        }
        mDS[i,2]=sDS[i,2]+sDS[k,2]
        mDS[i,1]=sDS[i,1]
        mDS[i,6]=x
        mDS[i,7]=N
      }
    mDS=as.data.frame(mDS)
    names(mDS)=c(names(sDS),'cov','n')
    mDS$scenario=scenario
    DS=rbind(DS,mDS)
    }
  }
}
DS=melt(DS,id.vars=c('maf','sites','n','cov','scenario'))
names(DS)[6:7]=c('Method','r2')
DSr=DS

DS$scenario=factor(DS$scenario,levels=c('LCWGS_low','LCWGS_med','LCWGS_hi'),labels=c('Low Diversity\nHigh LD\n','Med Diversity\nMed LD\n','Med Diversity\nLow LD\n'))
DS$cov=factor(DS$cov,levels=c(1,2,4),labels=c('1x','2x','4x'))
DS$method=factor(DS$Method,levels=c('angsd_r2','beagle_r2','stitch_r2'),labels=c('No Imputation','Beagle','Stitch'))

gg=ggplot(data=DS,mapping=aes(x=maf,y=r2,col=method))+facet_rep_grid(cov~scenario,scales="free_y",switch = "y")+geom_line(aes(linetype = as.factor(n)),size=1.1)

############Figure S20###########
#png('~/USS/aryn/LCWGS/Imputation_r2_byMAF_ryb.png',height=7,width=9,units="in",res=300,bg = "transparent")
print(gg+scale_color_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+scale_linetype_manual('n',values=c('dotted','dotdash','dashed','longdash','solid'))+theme(text = element_text(size=16),axis.text = element_text(size=12),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=14),legend.title=element_text(size=14),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 16),strip.background = element_blank(),panel.spacing = unit(2, "lines"),legend.key.width = unit(3,"line"))+ylab(expression(italic(r^{"2"})))+xlab('MAF'))#+coord_cartesian(ylim = c(0.2,1))
#dev.off()

#################################### GC by MAF ###################################
gts=NULL
for (scenario in c('LCWGS_low','LCWGS_med','LCWGS_hi')){
setwd(paste(c('~/USS/aryn/LCWGS/',scenario),collapse=''))

for (n in c(25,100,250,500,1000)){
    for (cov in c(1,2,4)){
        gt=read.table(paste(c('CalledGT_byMAF_N',n,'_',cov,'x.txt'),collapse=''),header=T)
        gt$N=rep(n,nrow(gt))
        gt$Cov=rep(cov,nrow(gt))
        gt$MAF=seq(0.05,1,0.05)
        mDS=matrix(nrow=10,ncol=11)
        for (i in 1:10){
            k=20-(i-1)
            for (j in c(1:6)){
                mDS[i,j]=weighted.mean(c(gt[i,j],gt[k,j]),c(gt[i,7],gt[k,7]))
            }
            mDS[i,7]=gt[i,7]+gt[k,7]
            mDS[i,8]=gt[i,14]
            mDS[i,9]=gt[i,13]
            mDS[i,10]=gt[i,12]
            mDS[i,11]=scenario
            }
            gts=rbind(gts,as.data.frame(mDS))
        }
    }
}

gts[,1:10]=apply(gts[,1:10],2,as.numeric)
names(gts)=c(names(gt)[1:7],'maf','cov','n','scenario')
gtsr=gts

gts=melt(gts,id.vars=c('sites','maf','cov','n','scenario'))
gts$Metric=str_split_fixed(gts$variable,"_",Inf)[,2]
gts$Method=str_split_fixed(gts$variable,"_",Inf)[,1]
gts$Method=factor(gts$Method,levels=c('Angsd','Beagle','Stitch'),labels=c('No Imputation','Beagle','Stitch'))
gts$scenario=factor(gts$scenario,levels=c('LCWGS_low','LCWGS_med','LCWGS_hi'),labels=c('Low Diversity\nHigh LD\n','Med Diversity\nMed LD\n','Med Diversity\nLow LD\n'))
gts$cov=factor(gts$cov,levels=c(1,2,4),labels=c('1x','2x','4x'))

############Figure S21###########
gg=ggplot(data=gts[which(gts$Metric=="corr"),],mapping=aes(x=maf,y=value,col=Method))+facet_rep_grid(cov~scenario,scales="free_y",switch = "y")+geom_line(aes(linetype = as.factor(n)),size=1.1)

#png('~/USS/aryn/LCWGS/ImputationGC_byMAF_ryb.png',height=7,width=9,units="in",res=300,bg = "transparent")
print(gg+scale_color_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+scale_linetype_manual('n',values=c('dotted','dotdash','dashed','longdash','solid'))+theme(text = element_text(size=16),axis.text = element_text(size=12),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=14),legend.title=element_text(size=14),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 16),strip.background = element_blank(),panel.spacing = unit(2, "lines"),legend.key.width = unit(3,"line"))+ylab('Genotype concordance')+xlab('MAF'))
#dev.off()

############Figure S22###########
gg=ggplot(data=gts[which(gts$Metric=="called"),],mapping=aes(x=maf,y=value,col=Method))+facet_rep_grid(cov~scenario,scales="free_y",switch = "y")+geom_line(aes(linetype = as.factor(n)),size=1.1)

#png('~/USS/aryn/LCWGS/ImputationCalled_byMAF_ryb.png',height=7,width=9,units="in",res=300,bg = "transparent")
print(gg+scale_color_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+scale_linetype_manual('n',values=c('dotted','dotdash','dashed','longdash','solid'))+theme(text = element_text(size=16),axis.text = element_text(size=12),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=14),legend.title=element_text(size=14),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 16),strip.background = element_blank(),panel.spacing = unit(2, "lines"),legend.key.width = unit(3,"line"))+ylab('Proportion Called')+xlab('MAF'))
#dev.off()


##################################### r2, called and GC by depth ##########################################
corrs=NULL
for (scenario in c('LCWGS_low','LCWGS_med','LCWGS_hi')){
setwd(paste(c('~/USS/aryn/LCWGS/',scenario),collapse=''))

gts=NULL
for (n in c(100,250,500,1000)){#25,
    for (cov in c(1,2,4)){
    gt=read.table(paste(c('CalledGT_byMAF_N',n,'_',cov,'x.txt'),collapse=''),header=T)
    gt$N=rep(n,nrow(gt))
    gt$Cov=rep(cov,nrow(gt))
    gt$MAF=seq(0.05,1,0.05)
    gts=rbind(gts,gt)
    }
}
gts$Ncov=apply(cbind(gts$N,gts$Cov),1,paste0,collapse='_')

#weighted mean of correct
gtss=gts[which(gts$MAF>0.1 & gts$MAF<0.9),]
Beagle=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Beagle_corr, x$sites*x$Beagle_called))[1:length(unique(gtss$Ncov))]
Stitch=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Stitch_corr, x$sites*x$Stitch_called))[1:length(unique(gtss$Ncov))]
Angsd=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Angsd_corr, x$sites*x$Angsd_called))[1:length(unique(gtss$Ncov))]
corr1=t(rbind(Beagle,Stitch,Angsd))
corr1=melt(corr1)
names(corr1)=c('Ncov','Method','Correct')
corr1$N=as.integer(sapply(as.character(corr1$Ncov),function(x) str_split_fixed(x,'_',Inf)[1]))
corr1$Cov=as.integer(sapply(as.character(corr1$Ncov),function(x) str_split_fixed(x,'_',Inf)[2]))

#weighted mean of called
gtss=gts[which(gts$MAF>0.1 & gts$MAF<0.9),]
Beagle=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Beagle_called, x$sites))[1:length(unique(gtss$Ncov))]
Stitch=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Stitch_called, x$sites))[1:length(unique(gtss$Ncov))]
Angsd=by(gtss, gtss$Ncov, function(x) weighted.mean(x$Angsd_called, x$sites))[1:length(unique(gtss$Ncov))]
corr2=t(rbind(Beagle,Stitch,Angsd))
corr2=melt(corr2)
names(corr2)=c('Ncov','Method','Called')
corr2$N=as.integer(sapply(as.character(corr2$Ncov),function(x) str_split_fixed(x,'_',Inf)[1]))
corr2$Cov=as.integer(sapply(as.character(corr2$Ncov),function(x) str_split_fixed(x,'_',Inf)[2]))

gts=NULL
for (n in c(100,250,500,1000)){#25,
    for (cov in c(1,2,4)){
    gt=read.table(paste(c('DosageR2_N',n,'_',cov,'x.txt'),collapse=''),header=T)
    gt$N=rep(n,nrow(gt))
    gt$Cov=rep(cov,nrow(gt))
    gts=rbind(gts,gt)
    }
}
gts$Ncov=apply(cbind(gts$N,gts$Cov),1,paste0,collapse='_')

gtss=gts[which(gts$maf>0.05 & gts$maf<0.95),]
Beagle=by(gtss, gtss$Ncov, function(x) weighted.mean(x$beagle_r2, x$sites))[1:length(unique(gtss$Ncov))]
Stitch=by(gtss, gtss$Ncov, function(x) weighted.mean(x$stitch_r2, x$sites))[1:length(unique(gtss$Ncov))]
Angsd=by(gtss, gtss$Ncov, function(x) weighted.mean(x$angsd_r2, x$sites))[1:length(unique(gtss$Ncov))]
corr3=t(rbind(Beagle,Stitch,Angsd))
corr3=melt(corr3)
names(corr3)=c('Ncov','Method','r2')
corr3$N=as.integer(sapply(as.character(corr3$Ncov),function(x) str_split_fixed(x,'_',Inf)[1]))
corr3$Cov=as.integer(sapply(as.character(corr3$Ncov),function(x) str_split_fixed(x,'_',Inf)[2]))

corr=merge(corr1,corr2,by=c('Ncov','Method','N','Cov'))
corr=merge(corr,corr3,by=c('Ncov','Method','N','Cov'))
corr$Scenario=scenario
corrs=rbind(corrs,corr)
}
corr=corrs

corrs=melt(corr,id.vars=c('Ncov','Method','N','Cov','Scenario'))
names(corrs)[6]='Metric'

corrs$Scenario=factor(corrs$Scenario,levels=c('LCWGS_low','LCWGS_med','LCWGS_hi'),labels=c('Low Diversity\nHigh LD\n','Med Diversity\nMed LD\n','Med Diversity\nLow LD\n'))
corrs$Metric=factor(corrs$Metric,levels=c('r2','Correct','Called'),labels=c(expression(italic(r^{"2"})),'GC','Called'))
corrs$Method=factor(corrs$Method,levels=c('Angsd','Beagle','Stitch'),labels=c('No Imputation','Beagle','Stitch'))
corrs$Cov[which(corrs$Cov==2)]=2.5

#############Figure S23##################
#png('~/USS/aryn/LCWGS/Imputation_allscenarios_ryb.png',height=7,width=9,units="in",res=300,bg = "transparent")

gg=ggplot(data=corrs,mapping=aes(x=Cov,y=value,col=Method))+facet_rep_grid(Metric~Scenario,scales="free_y",labeller = labeller(.rows=label_parsed),switch = "y")+geom_line(aes(linetype = as.factor(N)),size=1.1)

print(gg+scale_color_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+scale_linetype_manual('n',values=c('dotted','dotdash','longdash','solid'))+theme(text = element_text(size=20),axis.text = element_text(size=18),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=16),legend.title=element_text(size=18),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 18),strip.background = element_blank(),axis.title.y=element_blank(),panel.spacing = unit(2, "lines"),legend.key.width = unit(3,"line"))+ scale_x_continuous(breaks=c(1,2.5,4),labels=c(1,2,4))+xlab('Coverage'))

#dev.off()


#############Figure 9##################
#png('~/USS/aryn/LCWGS/Imputation_allscenarios_ryb_justDS.png',height=3.5,width=9,units="in",res=300,bg = "transparent")

gg=ggplot(data=corrs[which(corrs$Metric==levels(corrs$Metric)[1]),],mapping=aes(x=Cov,y=value,col=Method))+facet_rep_grid(.~Scenario,scales="free_y",labeller = labeller(.rows=label_parsed),switch = "y")+geom_line(aes(linetype = as.factor(N)),size=1.1)

print(gg+scale_color_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+scale_linetype_manual('n',values=c('dotted','dotdash','longdash','solid'))+theme(text = element_text(size=20),axis.text = element_text(size=14),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=14),legend.title=element_text(size=16),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 14),strip.background = element_blank(),panel.spacing = unit(2, "lines"),legend.key.width = unit(3,"line"))+ scale_x_continuous(breaks=c(1,2.5,4),labels=c(1,2,4))+xlab('Coverage')+ylab(expression(italic(r)^2)))

#dev.off()


############################## MAF improvement by imputation ##############################
maf=read.table('~/USS/aryn/LCWGS/MAFimputationImprovement.txt',header=T)
maf$Stitch_Diff=maf$Stitch_r2-maf$Angsd_r2
maf$Beagle_Diff=maf$Beagle_r2-maf$Angsd_r2

maf=melt(maf,id.vars=c('N','Cov','Scenario'))
maf$Method=str_split_fixed(maf$variable,'_',Inf)[,1]
maf$Metric=str_split_fixed(maf$variable,'_',Inf)[,2]
maf=maf[,c(1:3,5:7)]

maf$Scenario=factor(maf$Scenario,levels=c('LCWGS_low','LCWGS_med','LCWGS_hi'),labels=c('Low Diversity\nHigh LD\n','Med Diversity\nMed LD\n','Med Diversity\nLow LD\n'))
maf$Metric=factor(maf$Metric,levels=c('r2','RMSE','Diff'),labels=c(expression(italic(r^{"2"})),'RMSE',"Difference"))
maf$Cov=factor(maf$Cov,levels=c(1,2,4))
maf$Method=factor(maf$Method,levels=c('Angsd','Beagle','Stitch'))

#############Figure S24##################
#png('~/USS/aryn/LCWGS/Imputation_MAF_r2diff_ryb.png',height=10,width=9,units="in",res=300,bg = "transparent")
gg=ggplot(data=maf[which(maf$Metric=='Difference' & maf$N>25),],mapping=aes(x=Cov,y=value,fill=Method))+facet_rep_grid(N~Scenario,scales="free_y",switch = "y")+geom_bar(stat='identity',position='dodge')

gg+scale_fill_manual('Method',values=c('#941100','#FFB900','#005493'),breaks = c('Beagle','Stitch','No Imputation'))+xlab('Coverage')+theme(text = element_text(size=20),axis.text = element_text(size=18),panel.background = element_blank(),axis.line = element_line(colour = "black", size = 0.75, linetype = "solid"),axis.ticks=element_line(colour = "black", size = 0.75, linetype = "solid"),legend.position='right',legend.text=element_text(size=16),legend.title=element_text(size=18),legend.key = element_blank(),strip.placement = "outside",strip.text = element_text(size = 18),strip.background = element_blank(),panel.spacing = unit(2, "lines"))+ylab(expression(paste("Change in ",italic(r^{"2"})," with imputation")))+scale_y_continuous(labels = function(x) format(x, scientific = F))
#dev.off()
```
