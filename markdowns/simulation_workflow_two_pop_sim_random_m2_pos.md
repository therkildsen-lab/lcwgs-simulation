Simulation workflow for two-population simulation with random mutations under selection
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
    -   [Get SNP lists](#get-snp-lists)
    -   [Get shell script for per population MAF estimation](#get-shell-script-for-per-population-maf-estimation)
    -   [Run the shell script for MAF](#run-the-shell-script-for-maf)
    -   [Get Fst between the two populations](#get-fst-between-the-two-populations)
    -   [Run the shell script for Fst](#run-the-shell-script-for-fst)

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
if [ ! -d /workdir/lcwgs-simulation/two_pop_sim/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/two_pop_sim/rep_$REP_ID
  cd /workdir/lcwgs-simulation/two_pop_sim/rep_$REP_ID
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
  -d \"OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim.slim"
write_lines(shell_script, "../shell_scripts/two_pop_sim.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/two_pop_sim.sh $k > '/workdir/lcwgs-simulation/nohups/two_pop_sim_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_'$REP_ID'/'
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
write_lines(shell_script, "../shell_scripts/art_two_pop_sim.sh")
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_two_pop_sim.sh $k > '/workdir/lcwgs-simulation/nohups/art_two_pop_sim_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_'$REP_ID'/'
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
write_lines(shell_script, "../shell_scripts/merge_sort_subsample_two_pop_sim.sh")
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_two_pop_sim.sh $k > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_two_pop_sim_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
for (population in 1:2) {
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80,160)){
      for (i in 1:sample_size){
        if (i==1){
          write_lines(paste0("/workdir/lcwgs-simulation/two_pop_sim/rep_1/bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0("../two_pop_sim/rep_1/sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"))
        } else {
          write_lines(paste0("/workdir/lcwgs-simulation/two_pop_sim/rep_1/bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0("../two_pop_sim/rep_1/sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"), append = T)
        }
        if (i==1 & population==1){
          write_lines(paste0("/workdir/lcwgs-simulation/two_pop_sim/rep_1/bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0("../two_pop_sim/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"))
        } else {
          write_lines(paste0("/workdir/lcwgs-simulation/two_pop_sim/rep_1/bam/p", population, "_sample_", i, "_sorted_", coverage, "x.bam"), paste0("../two_pop_sim/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"), append = T)
        }
      }
    }
  }
}
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Get shell script for SNP calling
--------------------------------

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
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
    -setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 \\
    >& '/workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.log' &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/snp_calling_two_pop_sim.sh")
```

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_two_pop_sim.sh \
> /workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim.nohup &
```

Get SNP lists
-------------

``` bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
N_JOB_MAX=16
COUNT=0
## Create SNP lists
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    gunzip -c $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > \
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
    /workdir/programs/angsd0.931/angsd/angsd sites index \
    $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' &
    ## Submit sixteen jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done
```

Get shell script for per population MAF estimation
--------------------------------------------------

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
N_JOB_MAX=4
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
      -P 8 \\
      -setMinDepth 1 -minInd 1 -minQ 20 \\
      >& '/workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.log' & 
      ## Submit two jobs at a time
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_JOB_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done"
write_lines(shell_script, "../shell_scripts/get_maf_two_pop_sim.sh")
```

Run the shell script for MAF
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_maf_two_pop_sim.sh \
> /workdir/lcwgs-simulation/nohups/get_maf_two_pop_sim.nohup &
```

Get Fst between the two populations
-----------------------------------

``` r
shell_script <-"#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)                
        /workdir/programs/angsd0.931/angsd/misc/realSFS \\
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS'
        # Get alpha_beta
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst index \\
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
    -sfs $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' \\
    -fstout $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta'
        # Get alpha_beta.txt
        /workdir/programs/angsd0.931/angsd/misc/realSFS fst print \\
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.fst.idx' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt'
        # Get fst
        awk '{ print $0 \"\t\" $3 / $4 }' \\
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' \\
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.fst'
  done
done"
write_lines(shell_script, "../shell_scripts/get_fst_two_pop_sim.sh")
```

Run the shell script for Fst
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/get_fst_two_pop_sim.sh \
> /workdir/lcwgs-simulation/nohups/get_fst_two_pop_sim.nohup &
```
