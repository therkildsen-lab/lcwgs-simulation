Simulation workflow for two-population simulation with mutations under selection at fixed positions
================

-   [Two populations with divergent selection](#two-populations-with-divergent-selection)
    -   [Create a shell script to run SLiM with nohup](#create-a-shell-script-to-run-slim-with-nohup)
    -   [Run the shell script for SLiM simulation on server](#run-the-shell-script-for-slim-simulation-on-server)
    -   [Create a shell script to run ART with nohup](#create-a-shell-script-to-run-art-with-nohup)
    -   [Run the shell script for ART simulation on server](#run-the-shell-script-for-art-simulation-on-server)
    -   [Merge, sort, and subsample bam files](#merge-sort-and-subsample-bam-files)
    -   [Run the shell script for merging, sorting, and subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling)
    -   [Make bam lists](#make-bam-lists)
    -   [Index the ancestral fasta file](#index-the-ancestral-fasta-file)
    -   [Get shell script for SNP calling](#get-shell-script-for-snp-calling)
    -   [Run the shell script for SNP calling](#run-the-shell-script-for-snp-calling)
    -   [Get shell script writing SNP lists](#get-shell-script-writing-snp-lists)
    -   [Run shell script to get SNP list](#run-shell-script-to-get-snp-list)
    -   [Get shell script for per population MAF estimation](#get-shell-script-for-per-population-maf-estimation)
    -   [Run the shell script for MAF](#run-the-shell-script-for-maf)
    -   [Get Fst between the two populations](#get-fst-between-the-two-populations)
    -   [Run the shell script for Fst](#run-the-shell-script-for-fst)
-   [Two populations with divergent selection, with lower selection](#two-populations-with-divergent-selection-with-lower-selection)
    -   [Create a shell script to run SLiM with nohup](#create-a-shell-script-to-run-slim-with-nohup-1)
    -   [Run the shell script for SLiM simulation on server](#run-the-shell-script-for-slim-simulation-on-server-1)
    -   [Run the shell script for ART simulation on server](#run-the-shell-script-for-art-simulation-on-server-1)
    -   [Run the shell script for merging, sorting, and subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling-1)
    -   [Make bam lists](#make-bam-lists-1)
    -   [Index the ancestral fasta file](#index-the-ancestral-fasta-file-1)
    -   [Run the shell script for SNP calling](#run-the-shell-script-for-snp-calling-1)
    -   [Run the shell script to get SNP lists](#run-the-shell-script-to-get-snp-lists)
    -   [Run the shell script for MAF](#run-the-shell-script-for-maf-1)
    -   [Run the shell script for Fst](#run-the-shell-script-for-fst-1)
-   [Two populations with divergent selection, with lower selection and lower recombination](#two-populations-with-divergent-selection-with-lower-selection-and-lower-recombination)
    -   [Create a shell script to run SLiM with nohup](#create-a-shell-script-to-run-slim-with-nohup-2)
    -   [Run the shell script for SLiM simulation on server](#run-the-shell-script-for-slim-simulation-on-server-2)
    -   [Run the shell script for ART simulation on server](#run-the-shell-script-for-art-simulation-on-server-2)
    -   [Run the shell script for merging, sorting, and subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling-2)
    -   [Make bam lists](#make-bam-lists-2)
    -   [Index the ancestral fasta file](#index-the-ancestral-fasta-file-2)
    -   [Run the shell script for SNP calling](#run-the-shell-script-for-snp-calling-2)
    -   [Run the shell script to get SNP lists](#run-the-shell-script-to-get-snp-lists-1)
    -   [Run the shell script for MAF](#run-the-shell-script-for-maf-2)
    -   [Run the shell script for Fst](#run-the-shell-script-for-fst-2)
-   [Two populations with divergent selection, with lower selection and lower migration rate](#two-populations-with-divergent-selection-with-lower-selection-and-lower-migration-rate)
    -   [Create a shell script to run SLiM with nohup](#create-a-shell-script-to-run-slim-with-nohup-3)
    -   [Run the shell script for SLiM simulation on server](#run-the-shell-script-for-slim-simulation-on-server-3)
    -   [Run the shell script for ART simulation on server](#run-the-shell-script-for-art-simulation-on-server-3)
    -   [Run the shell script for merging, sorting, and subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling-3)
    -   [Make bam lists](#make-bam-lists-3)
    -   [Index the ancestral fasta file](#index-the-ancestral-fasta-file-3)
    -   [Run the shell script for SNP calling](#run-the-shell-script-for-snp-calling-3)
    -   [Run the shell script to get SNP lists](#run-the-shell-script-to-get-snp-lists-2)
    -   [Run the shell script for MAF](#run-the-shell-script-for-maf-3)
    -   [Run the shell script for Fst](#run-the-shell-script-for-fst-3)

``` r
library(tidyverse)
```

Two populations with divergent selection
========================================

Create a shell script to run SLiM with nohup
--------------------------------------------

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
  -d MUTATION_RATE=2e-7 \\
  -d SELECTION_COEFF=0.25 \\
  -d MIGRATION_RATE=0.05 \\
  -d REC_RATE=1e-7 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=1000 \\
  -d N_M2=6 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_fixed_m2_pos.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

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

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

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

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_'$k'.nohup' &
done
```

Make bam lists
--------------

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

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Get shell script for SNP calling
--------------------------------

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=6
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    ## SNP calling
    /workdir/programs/angsd0.931/angsd/angsd \\
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
    -anc $BASE_DIR'slim/ancestral.fasta' \\
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
    -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
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

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos.nohup &
```

Get shell script writing SNP lists
----------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=16
COUNT=0
## Create SNP lists
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    gunzip -c $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > \\
    $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
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
    $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
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

Run shell script to get SNP list
--------------------------------

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'
```

Get shell script for per population MAF estimation
--------------------------------------------------

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=16
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
      for POP in {1,2}; do
      ## Get MAF
      /workdir/programs/angsd0.931/angsd/angsd \\
      -b $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -anc $BASE_DIR'slim/ancestral.fasta' \\
      -out $BASE_DIR'angsd/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
      -sites $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
      -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
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

Run the shell script for MAF
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos.nohup &
```

Get Fst between the two populations
-----------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
# Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)                
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS \\
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' &
  done
done
wait
# Get alpha_beta
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst index \\
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    -sfs $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' \\
    -fstout $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta' &
  done
done
wait
# Get alpha_beta.txt
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst print \\
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.fst.idx' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' &
  done
done
wait
# Get fst
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
        awk '{ print $0 \"\t\" $3 / $4 }' \\
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.fst' &
  done
done"
write_lines(shell_script, "../shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh")
```

Run the shell script for Fst
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos.nohup &
```

Two populations with divergent selection, with lower selection
==============================================================

Create a shell script to run SLiM with nohup
--------------------------------------------

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_$REP_ID
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
  -d MUTATION_RATE=2e-7 \\
  -d SELECTION_COEFF=0.10 \\
  -d MIGRATION_RATE=0.05 \\
  -d REC_RATE=1e-7 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=1000 \\
  -d N_M2=6 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_fixed_m2_pos.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_lower_s.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_lower_s.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_lower_s_'$k'.nohup' &
done
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_lower_s_'$k'.nohup' &
done
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_lower_s_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
make_bam_lists("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_1/")
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos_lower_s.nohup &
```

Run the shell script to get SNP lists
-------------------------------------

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/'
```

Run the shell script for MAF
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos_lower_s.nohup &
```

Run the shell script for Fst
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos_lower_s.nohup &
```

Two populations with divergent selection, with lower selection and lower recombination
======================================================================================

Create a shell script to run SLiM with nohup
--------------------------------------------

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
  -d MUTATION_RATE=2e-7 \\
  -d SELECTION_COEFF=0.10 \\
  -d MIGRATION_RATE=0.05 \\
  -d REC_RATE=2e-8 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=1000 \\
  -d N_M2=6 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_fixed_m2_pos.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_r.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_lower_s_lower_r_'$k'.nohup' &
done
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_lower_s_lower_r_'$k'.nohup' &
done
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_lower_s_lower_r_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
make_bam_lists("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

Run the shell script to get SNP lists
-------------------------------------

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'
```

Run the shell script for MAF
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

Run the shell script for Fst
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos_lower_s_lower_r.nohup &
```

Two populations with divergent selection, with lower selection and lower migration rate
=======================================================================================

Create a shell script to run SLiM with nohup
--------------------------------------------

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/rep_$REP_ID
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
  -d MUTATION_RATE=2e-7 \\
  -d SELECTION_COEFF=0.05 \\
  -d MIGRATION_RATE=0.01 \\
  -d REC_RATE=1e-7 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=1000 \\
  -d N_M2=6 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_fixed_m2_pos.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_m.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim_fixed_m2_pos_lower_s_lower_m.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_fixed_m2_pos_lower_s_lower_m_'$k'.nohup' &
done
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/' > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_fixed_m2_pos_lower_s_lower_m_'$k'.nohup' &
done
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim_fixed_m2_pos.sh $k '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_fixed_m2_pos_lower_s_lower_m_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
make_bam_lists("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/rep_1/")
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_fixed_m2_pos_lower_s_lower_m.nohup &
```

Run the shell script to get SNP lists
-------------------------------------

``` bash
bash /workdir/lcwgs-simulation/shell_scripts/get_snp_list_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/'
```

Run the shell script for MAF
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/' \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim_fixed_m2_pos_lower_s_lower_m.nohup &
```

Run the shell script for Fst
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim_fixed_m2_pos.sh 1 '/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_m/' \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim_fixed_m2_pos_lower_s_lower_m.nohup &
```
