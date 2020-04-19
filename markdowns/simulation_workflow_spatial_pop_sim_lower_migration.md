Simulation workflow for two dimensional spatial populations with lower migration rate
================

-   [Two dimensional spatial populations with lower migration rate](#two-dimensional-spatial-populations-with-lower-migration-rate)
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
    -   [Get shell script for PCA](#get-shell-script-for-pca)
    -   [Run the shell script for PCA](#run-the-shell-script-for-pca)

``` r
library(tidyverse)
```

Two dimensional spatial populations with lower migration rate
=============================================================

Create a shell script to run SLiM with nohup
--------------------------------------------

The only parameter that I've changed is the migration rate.

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/spatial_pop_sim_lower_m/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/spatial_pop_sim_lower_m/rep_$REP_ID
  cd /workdir/lcwgs-simulation/spatial_pop_sim_lower_m/rep_$REP_ID
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
  -d MUTATION_RATE=20e-8 \\
  -d REC_RATE=50e-8 \\
  -d CHR_LENGTH=30000000 \\
  -d META_POP_SIDE=3 \\
  -d POP_SIZE=500 \\
  -d MIGRATION_RATE=0.0005 \\
  -d SAMPLE_SIZE=160 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/spatial_pop_sim.slim"
write_lines(shell_script, "../shell_scripts/spatial_pop_sim_lower_m.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/spatial_pop_sim_lower_m.sh $k > '/workdir/lcwgs-simulation/nohups/spatial_pop_sim_lower_m_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=9
## Generate sam files
COUNT=0
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      /workdir/programs/art_bin_MountRainier/art_illumina \\
      -ss HS25 \\
      -sam \\
      -i $OUT_DIR'fasta/p'$j'_sample'$i'_genome'$k'.fasta' \\
      -p \\
      -na \\
      -l 150 \\
      -f 10 \\
      -m 500 \\
      -s 75 \\
      -rs $(($REP_ID*10000+$i)) \\
      -o $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k &
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
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      samtools view -bS -F 4 $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.sam' \\
      > $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.bam' &
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
for i in {1..80}; do
  for j in {1..9}; do
    for k in {1..2}; do
      rm $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'.sam'
      mv $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'1.fq' $OUT_DIR'fastq/'
      mv $OUT_DIR'bam/p'$j'_sample'$i'_genome'$k'2.fq' $OUT_DIR'fastq/'
    done
  done
done"
write_lines(shell_script, "../shell_scripts/art_spatial_pop_sim.sh")
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_spatial_pop_sim.sh $k '/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/' > '/workdir/lcwgs-simulation/nohups/art_spatial_pop_sim_lower_m_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
OUT_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_CORE_MAX=9
## merge
COUNT=0
for k in {1..80}; do
  for i in {1..9}; do
    samtools merge $OUT_DIR'bam/p'$i'_sample'$k'.bam' \\
    $OUT_DIR'bam/p'$i'_sample'$k'_genome1.bam' \\
    $OUT_DIR'bam/p'$i'_sample'$k'_genome2.bam' &
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
for k in {1..80}; do
  for i in {1..9}; do
    samtools sort -o $OUT_DIR'bam/p'$i'_sample'$k'_sorted.bam' $OUT_DIR'bam/p'$i'_sample'$k'.bam' &
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
for j in {0.125,0.25,0.5,1,2,4}; do
  for k in {1..80}; do
    for i in {1..9}; do
      samtools view \\
      -s `awk -v j=$j 'BEGIN { print j / 20 }'` \\
      -b $OUT_DIR'bam/p'$i'_sample'$k'_sorted.bam' \\
      > $OUT_DIR'bam/p'$i'_sample'$k'_sorted_'$j'x.bam' &
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
for k in {1..80}; do
  for i in {1..9}; do
    #rm $OUT_DIR'bam/p'$i'_sample'$k'_genome1.bam' 
    #rm $OUT_DIR'bam/p'$i'_sample'$k'_genome2.bam' 
    rm $OUT_DIR'bam/p'$i'_sample'$k'.bam'
  done
done"
write_lines(shell_script, "../shell_scripts/merge_sort_subsample_spatial_pop_sim.sh")
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_spatial_pop_sim.sh $k '/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/' > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_spatial_pop_sim_lower_m_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
make_bam_lists <- function(basedir){
  i=1
  for (population in 1:9) {
    for (coverage in c(0.125,0.25,0.5,1,2,4)){
      for (sample_size in c(5,10,20,40,80)){
        for (i in 1:sample_size){
          if (i==1){
            write_lines(paste0(basedir, "bam/p", population, "_sample", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"))
          } else {
            write_lines(paste0(basedir, "bam/p", population, "_sample", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_p", population, "_", sample_size, "_", coverage, "x.txt"), append = T)
          }
          if (i==1 & population==1){
            write_lines(paste0(basedir, "bam/p", population, "_sample", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"))
          } else {
            write_lines(paste0(basedir, "bam/p", population, "_sample", i, "_sorted_", coverage, "x.bam"), paste0(basedir, "sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"), append = T)
          }
        }
      }
    }
  }
}
make_bam_lists("/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/rep_1/")
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Get shell script for SNP calling
--------------------------------

After some experimentation, it turns out that a MAF filter is very important to have.

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
N_JOB_MAX=5
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80}; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
    ## SNP calling
    /workdir/programs/angsd0.931/angsd/angsd \\
    -b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
    -anc $BASE_DIR'slim/ancestral.fasta' \\
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
    -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
    -doCounts 1 -doDepth 1 -dumpCounts 1 \\
    -doIBS 2 -makematrix 1 -doCov 1 \\
    -P 2 -SNP_pval 1e-6 -rmTriallelic 1e-6 \\
    -setMinDepth 2 -minInd 1 -minMaf 0.05 -minQ 20 &
    ## Submit two jobs at a time
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_JOB_MAX ]; then
     wait
     COUNT=0
    fi
  done
done"
write_lines(shell_script, "../shell_scripts/snp_calling_spatial_pop_sim.sh")
```

Run the shell script for SNP calling
------------------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_spatial_pop_sim.sh 1 '/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/' \
> /workdir/lcwgs-simulation/nohups/snp_calling_spatial_pop_sim_lower_m.nohup &
```

Get shell script for PCA
------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR_BASE=$2
BASE_DIR=$OUT_DIR_BASE'rep_'$REP_ID'/'
for SAMPLE_SIZE in {5,10,20,40,80}; do
  for COVERAGE in {0.125,0.25,0.5,1,2,4}; do
      ## Run PCAngsd
      python /workdir/programs/pcangsd/pcangsd.py \\
      -beagle $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.beagle.gz' \\
      -minMaf 0.05 \\
      -threads 16 \\
      -o $BASE_DIR'angsd/pcagnsd_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x'
  done
done"
write_lines(shell_script, "../shell_scripts/run_pcangd_spatial_pop_sim.sh")
```

Run the shell script for PCA
----------------------------

``` bash
nohup bash /workdir/lcwgs-simulation/shell_scripts/run_pcangd_spatial_pop_sim.sh 1 '/workdir/lcwgs-simulation/spatial_pop_sim_lower_m/' \
> /workdir/lcwgs-simulation/nohups/run_pcangd_spatial_pop_sim_lower_m.nohup &
```
