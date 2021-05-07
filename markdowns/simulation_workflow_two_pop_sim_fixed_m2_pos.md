Simulation workflow for two-population simulation with mutations under
selection at fixed positions
================

  - [Two populations with divergent selection ( Ne\~50,000 in each
    population)](#two-populations-with-divergent-selection-ne50000-in-each-population)
      - [Create a shell script to run the burnin step with
        SLiM](#create-a-shell-script-to-run-the-burnin-step-with-slim)
      - [Run the shell script for the burnin step with
        SLiM](#run-the-shell-script-for-the-burnin-step-with-slim)
      - [Create a shell script to run the selection step with
        SLiM](#create-a-shell-script-to-run-the-selection-step-with-slim)
      - [Run the shell script for the selection step with
        SLiM](#run-the-shell-script-for-the-selection-step-with-slim)
      - [Create a shell script to run ART with
        nohup](#create-a-shell-script-to-run-art-with-nohup)
      - [Run the shell script for ART simulation on
        server](#run-the-shell-script-for-art-simulation-on-server)
      - [Merge, sort, and subsample bam
        files](#merge-sort-and-subsample-bam-files)
      - [Run the shell script for merging, sorting, and
        subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling)
      - [Make bam lists](#make-bam-lists)
      - [Index the ancestral fasta
        file](#index-the-ancestral-fasta-file)
      - [Get Fst with Samtool’s GL
        model](#get-fst-with-samtools-gl-model)
      - [Theta and neutrality stats with Samtools GL model in population
        1](#theta-and-neutrality-stats-with-samtools-gl-model-in-population-1)
      - [Testing GATK’s GL model](#testing-gatks-gl-model)
      - [Get shell script for running selection scan using
        PCAngsd](#get-shell-script-for-running-selection-scan-using-pcangsd)
      - [Run the shell script for selection scan using
        PCAngsd](#run-the-shell-script-for-selection-scan-using-pcangsd)
  - [Two populations with divergent selection, with smaller population
    size ( Ne\~10,000 in each
    population)](#two-populations-with-divergent-selection-with-smaller-population-size-ne10000-in-each-population)
      - [Create a shell script to run the burnin step with
        SLiM](#create-a-shell-script-to-run-the-burnin-step-with-slim-1)
      - [Run the shell script for the burnin step with
        SLiM](#run-the-shell-script-for-the-burnin-step-with-slim-1)
      - [Create a shell script to run the selection step with
        SLiM](#create-a-shell-script-to-run-the-selection-step-with-slim-1)
      - [Run the shell script for the selection step with
        SLiM](#run-the-shell-script-for-the-selection-step-with-slim-1)
      - [Run the shell script for ART simulation on
        server](#run-the-shell-script-for-art-simulation-on-server-1)
      - [Run the shell script for merging, sorting, and
        subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling-1)
      - [Make bam lists](#make-bam-lists-1)
      - [Index the ancestral fasta
        file](#index-the-ancestral-fasta-file-1)
      - [Run the shell script for SNP
        calling](#run-the-shell-script-for-snp-calling-2)
      - [Run the shell script to get SNP
        lists](#run-the-shell-script-to-get-snp-lists)
      - [Run the shell script for MAF](#run-the-shell-script-for-maf-2)
      - [Run the shell script for Fst](#run-the-shell-script-for-fst-2)
      - [Run the shell script for selection scan using
        PCAngsd](#run-the-shell-script-for-selection-scan-using-pcangsd-1)

``` r
library(tidyverse)
```

# Two populations with divergent selection ( Ne\~50,000 in each population)

## Create a shell script to run the burnin step with SLiM

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Run SLiM 
/workdir/programs/SLiM_build/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=100e-8 \\
  -d REC_RATE=250e-8 \\
  -d MIGRATION_RATE=0.01 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=500 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_burnin.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_burnin.sh")
```

## Run the shell script for the burnin step with SLiM

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_burnin.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_burnin_'$k'.nohup' &
done
```

## Create a shell script to run the selection step with SLiM

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Run SLiM 
/workdir/programs/SLiM_build/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=100e-9 \\
  -d REC_RATE=250e-9 \\
  -d MIGRATION_RATE=0.001 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=5000 \\
  -d SAMPLE_SIZE=1000 \\
  -d SELECTION_COEFF=0.08 \\
  -d N_M2=11 \\
  -d M2_FREQUENCY=1 \\
  -d \"ANCESTRAL_FASTA='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/ancestral.fasta'\" \\
  -d \"BURNIN_FILE='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/burnin_full_output.txt'\" \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_selection.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_selection.sh")
```

## Run the shell script for the selection step with SLiM

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_selection.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_selection_'$k'.nohup' &
done
```

## Create a shell script to run ART with nohup

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=28
## Generate sam files
COUNT=0
for i in {1..160}; do
  for j in {1..2}; do
    for k in {1..2}; do
      /workdir/programs/art_bin_MountRainier/art_illumina \\
      -ss HS25 \\
      -sam \\
      -i $OUT_DIR'fasta/p'$j'_derived_'$i'_'$k'.fasta' \\
      -p \\
      -na \\
      -l 150 \\
      -f 10 \\
      -m 500 \\
      -s 75 \\
      -rs $(($REP_ID*10000+$i)) \\
      -o $OUT_DIR'bam/p'$j'_derived_'$i'_'$k &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
        wait
        COUNT=0
      fi
    done
  done
done
wait
## Generate bam files
COUNT=0
for i in {1..160}; do
  for j in {1..2}; do
    for k in {1..2}; do
      samtools view -bS -F 4 $OUT_DIR'bam/p'$j'_derived_'$i'_'$k'.sam' \\
      > $OUT_DIR'bam/p'$j'_derived_'$i'_'$k'.bam' &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
        wait
        COUNT=0
      fi
    done
  done
done
wait
## Sort bam files
for i in {1..320}; do
  for j in {1..2}; do
    for k in {1..2}; do
      rm $OUT_DIR'bam/p'$j'_derived_'$i'_'$k'.sam'
      mv $OUT_DIR'bam/p'$j'_derived_'$i'_'$k'1.fq' $OUT_DIR'fastq/'
      mv $OUT_DIR'bam/p'$j'_derived_'$i'_'$k'2.fq' $OUT_DIR'fastq/'
    done
  done
done"
write_lines(shell_script, "../shell_scripts/art_two_pop_sim_fixed_m2_pos.sh")
```

## Run the shell script for ART simulation on server

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_'$k'.nohup' &
done
```

## Merge, sort, and subsample bam files

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=40
## merge
COUNT=0
for k in {1..160}; do
  for i in {1..2}; do
    samtools merge $OUT_DIR'bam/p'$i'_sample_'$k'.bam' \\
    $OUT_DIR'bam/p'$i'_derived_'$k'_1.bam' \\
    $OUT_DIR'bam/p'$i'_derived_'$k'_2.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## sort
COUNT=0
for k in {1..160}; do
  for i in {1..2}; do
    samtools sort -o $OUT_DIR'bam/p'$i'_sample_'$k'_sorted.bam' $OUT_DIR'bam/p'$i'_sample_'$k'.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## subsample
COUNT=0
for j in {0.25,0.5,1,2,4,8}; do
  for k in {1..160}; do
    for i in {1..2}; do
      samtools view \\
      -s `awk -v j=$j 'BEGIN { print j / 20 }'` \\
      -b $OUT_DIR'bam/p'$i'_sample_'$k'_sorted.bam' \\
      > $OUT_DIR'bam/p'$i'_sample_'$k'_sorted_'$j'x.bam' &
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_CORE_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done
wait
## delete intermediate files
for k in {1..160}; do
  for i in {1..2}; do
    l=$(($k+160))
    #rm $OUT_DIR'bam/p'$i'_derived_'$k'.bam' 
    #rm $OUT_DIR'bam/p'$i'_derived_'$l'.bam' 
    rm $OUT_DIR'bam/p'$i'_sample_'$k'.bam'
  done
done"
write_lines(shell_script, "../shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh")
```

## Run the shell script for merging, sorting, and subsampling

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_'$k'.nohup' &
done
```

## Make bam lists

``` r
make_bam_lists <- function(basedir){
  i=1
  for (population in 1:2) {
    for (coverage in c(0.25,0.5,1,2,4,8)){
      for (sample_size in c(5,10,20,40,80,160)){
        for (i in 1:sample_size){
          if (i==1){
            write_lines(paste0(basedir, "bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"))
          } else {
            write_lines(paste0(basedir, "bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"), append = T)
          }
          if (i==1 & population==1){
            write_lines(paste0(basedir, "bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"))
          } else {
            write_lines(paste0(basedir, "bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"), append = T)
          }
        }
      }
    }
  }
}
make_bam_lists("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/")
```

## Index the ancestral fasta file

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

## Get Fst with Samtool’s GL model

#### Get shell script for SNP calling

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
GL_MODEL=${4:-1}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=6
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    ## SNP calling
    /workdir/programs/angsd0.931/angsd/angsd \\
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
    -anc $BASE_DIR'slim/ancestral.fasta' \\
    -out $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
    -dosaf 1 -GL $GL_MODEL -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
    -doCounts 1 -doDepth 1 -dumpCounts 1 \\
    -P 6 -SNP_pval 1e-6 -rmTriallelic 1e-6 \\
    -setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh")
```

#### Run the shell script for SNP calling

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos.nohup &
```

#### Get shell script writing SNP lists

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=16
COUNT=0
## Create SNP lists
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    gunzip -c $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > \\
    $BASE_DIR$SUBDIR'/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
    ## Submit sixteen jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
COUNT=0
## Index the SNP lists
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd sites index \\
    $BASE_DIR$SUBDIR'/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
    ## Submit sixteen jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh")
```

#### Run shell script to get SNP list

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'
```

#### Get shell script for per population MAF estimation

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
GL_MODEL=${4:-1}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=12
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
      for POP in {1,2}; do
      ## Get MAF
      /workdir/programs/angsd0.931/angsd/angsd \\
      -b $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -anc $BASE_DIR'slim/ancestral.fasta' \\
      -out $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
      -sites $BASE_DIR$SUBDIR'/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -dosaf 1 -GL $GL_MODEL -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
      -doCounts 1 -doDepth 1 -dumpCounts 1 \\
      -P 2 \\
      -setMinDepth 1 -minInd 1 -minQ 20 & 
      ## Submit two jobs at a time
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_JOB_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done"
write_lines(shell_script, "../shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh")
```

#### Run the shell script for MAF

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos.nohup &
```

#### Get Fst between the two populations

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
SUBDIR=${3:-angsd}
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=2
COUNT=0
# Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)                
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS \\
    $BASE_DIR$SUBDIR'/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR$SUBDIR'/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    -P 8 \\
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' & 
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
# Get alpha_beta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst index \\
    $BASE_DIR$SUBDIR'/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR$SUBDIR'/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    -P 8 \\
    -sfs $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' \\
    -fstout $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta' &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
wait
# Get alpha_beta.txt
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst print \\
    $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.fst.idx' \\
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' &
  done
done
wait
# Get fst
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        awk '{ print $0 \"\t\" $3 / $4 }' \\
    $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' \\
    > $BASE_DIR$SUBDIR'/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.fst' &
  done
done"
write_lines(shell_script, "../shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh")
```

#### Run the shell script for Fst

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos.nohup &
```

## Theta and neutrality stats with Samtools GL model in population 1

#### Create shell scripts to estimated theta and neutrality stats

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/}
SUBDIR=${3:-angsd}
GL_MODEL=${4:-1}
POP=${5:-1}

BASE_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=6

## Get saf file
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \\
      -bam $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -out $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites' \\
      -doSaf 1 \\
      -anc $BASE_DIR'slim/ancestral.fasta' \\
      -GL $GL_MODEL \\
      -P 6 \\
      -doCounts 1 \\
      -setMinDepth `awk \"BEGIN {print $SAMPLE_SIZE*$COVERAGE}\"` &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Get SFS from saf (this need to be broken up into two separate runs due to memory limiations)
COUNT=0
N_CORE_MAX=1
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/realSFS \\
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.saf.idx' \\
      -P 30 \\
      -tole 1e-08 \\
      -maxIter 1000 \\
      > $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.sfs' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Estimate theta
COUNT=0
N_CORE_MAX=6
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/angsd \\
      -bam $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -out $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites' \\
      -doThetas 1 \\
      -doSaf 1 \\
      -pest $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.sfs' \\
      -anc $BASE_DIR'slim/ancestral.fasta' \\
      -GL $GL_MODEL \\
      -P 6 \\
      -doCounts 1 \\
      -setMinDepth `awk \"BEGIN {print $SAMPLE_SIZE*$COVERAGE}\"` &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Print per-SNP theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat print \\
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \\
      > $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.tsv' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Do fixed window theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \\
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \\
      -win 10000 -step 10000 \\
      -outnames $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.windowed_thetas.idx' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
## Do per-chromosome average theta
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \\
      $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.thetas.idx' \\
      -outnames $BASE_DIR$SUBDIR'/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x_all_sites.average_thetas.idx' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/estimate_theta_two_populations.sh")
```

#### Run shell scripts to estimated theta and neutrality stats with Samtools’ GL model in p1

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/estimate_theta_two_populations.sh \
1 \
/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/ \
angsd \
1 \
1 \
> /workdir/lcwgs-simulation/nohups/estimate_theta_two_pop_sim_fixed_m2_pos_p1.nohup &
```

## Testing GATK’s GL model

For this step, I used GATK’s GL model (`-GL 2`).

#### Run the shell script for SNP calling

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
angsd_gatk \
2 \
> /workdir/lcwgs-simulation/nohups/snp_calling_gatk_two_pop_sim_fixed_m2_pos.nohup &
```

#### Run shell script to get SNP list

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
angsd_gatk
```

#### Run the shell script for MAF

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
angsd_gatk \
2 \
> /workdir/lcwgs-simulation/nohups/get_maf_gatk_two_pop_sim_fixed_m2_pos.nohup &
```

#### Run the shell script for Fst

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
angsd_gatk \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos.nohup &
```

## Get shell script for running selection scan using PCAngsd

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=3
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    ## Run PCA
    python /workdir/programs/pcangsd/pcangsd.py \\
    -beagle $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.beagle.gz' \\
    -selection \\
    -sites_save \\
    -minMaf 0.05 \\
    -threads 10 \\
    -iter 200 \\
    -maf_iter 200 \\
    -o $BASE_DIR'angsd/pcagnsd_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' &
    ## Submit three jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/run_pcangsd_two_pop_sim_fixed_m2_pos.sh")
```

## Run the shell script for selection scan using PCAngsd

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/run_pcangsd_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/run_pcangsd_two_pop_sim_fixed_m2_pos.nohup &
```

# Two populations with divergent selection, with smaller population size ( Ne\~10,000 in each population)

The same population size is simulated, but I’ve scaled down mutation
rate, recombination rate, migration rate. The selection coefficient is
unchanged. (Ignore the directory name. I’m just using it for
convenience.)

## Create a shell script to run the burnin step with SLiM

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Run SLiM 
/workdir/programs/SLiM_build/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=200e-9 \\
  -d REC_RATE=500e-9 \\
  -d MIGRATION_RATE=0.005 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=500 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_burnin.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r_burnin.sh")
```

## Run the shell script for the burnin step with SLiM

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r_burnin.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_lower_s_lower_r_burnin_'$k'.nohup' &
done
```

## Create a shell script to run the selection step with SLiM

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Run SLiM 
/workdir/programs/SLiM_build/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=200e-10 \\
  -d REC_RATE=500e-10 \\
  -d MIGRATION_RATE=0.0005 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=5000 \\
  -d SAMPLE_SIZE=1000 \\
  -d SELECTION_COEFF=0.08 \\
  -d N_M2=11 \\
  -d M2_FREQUENCY=1 \\
  -d \"ANCESTRAL_FASTA='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/slim/ancestral.fasta'\" \\
  -d \"BURNIN_FILE='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/slim/burnin_full_output.txt'\" \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_selection.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection.sh")
```

## Run the shell script for the selection step with SLiM

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection_'$k'.nohup' &
done
```

## Run the shell script for ART simulation on server

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_lower_s_lower_r_'$k'.nohup' &
done
```

## Run the shell script for merging, sorting, and subsampling

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_lower_s_lower_r_'$k'.nohup' &
done
```

## Make bam lists

``` r
make_bam_lists("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
```

## Index the ancestral fasta file

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

## Run the shell script for SNP calling

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

## Run the shell script to get SNP lists

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'
```

## Run the shell script for MAF

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

## Run the shell script for Fst

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

## Run the shell script for selection scan using PCAngsd

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/run_pcangsd_two_pop_sim_fixed_m2_pos.sh \
1 \
'/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/run_pcangsd_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```
