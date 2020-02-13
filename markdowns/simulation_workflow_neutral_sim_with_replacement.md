Simulation workflow for neutral simulation, sampling with replacement
================

-   [Neutral simulation, sampling with replacement](#neutral-simulation-sampling-with-replacement)
    -   [Copy fasta files from the original neutral simulation directory](#copy-fasta-files-from-the-original-neutral-simulation-directory)
    -   [Run the shell script for sampling with replacement](#run-the-shell-script-for-sampling-with-replacement)
    -   [Create a shell script to run ART with nohup](#create-a-shell-script-to-run-art-with-nohup)
    -   [Run the shell script for ART simulation on server](#run-the-shell-script-for-art-simulation-on-server)
    -   [Merge, sort, and subsample bam files](#merge-sort-and-subsample-bam-files)
    -   [Run the shell script for merging, sorting, and subsampling](#run-the-shell-script-for-merging-sorting-and-subsampling)
    -   [Make bam lists](#make-bam-lists)

``` r
library(tidyverse)
```

Neutral simulation, sampling with replacement
=============================================

Copy fasta files from the original neutral simulation directory
---------------------------------------------------------------

In this step, randomly sample these fasta files with replacement.

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID
  cd /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Copy fasta files
for I in {1..2000}; do
  INDEX=`shuf -i 1-2000 -n 1`
  echo $I'_'$INDEX
  cp '/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/fasta/derived_'$INDEX'.fasta' '/workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_'$REP_ID'/fasta/derived_'$I'.fasta'
done"
write_lines(shell_script, "../shell_scripts/sample_with_replacement.sh")
```

Run the shell script for sampling with replacement
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/sample_with_replacement.sh $k > '/workdir/lcwgs-simulation/nohups/sample_with_replacement_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/neutral_sim/}

OUT_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=40
## Generate sam files
COUNT=0
for i in {1..2000}; do
  /workdir/programs/art_bin_MountRainier/art_illumina \\
  -ss HS25 \\
  -sam \\
  -i $OUT_DIR'fasta/derived_'$i'.fasta' \\
  -p \\
  -na \\
  -l 150 \\
  -f 10 \\
  -m 500 \\
  -s 75 \\
  -rs $(($REP_ID*10000+$i)) \\
  -o $OUT_DIR'bam/derived_'$i &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## Generate bam files
COUNT=0
for i in {1..2000}; do
  samtools view -bS -F 4 $OUT_DIR'bam/derived_'$i'.sam' > $OUT_DIR'bam/derived_'$i'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## Sort bam files
for i in {1..2000}; do
  rm $OUT_DIR'bam/derived_'$i'.sam'
  mv $OUT_DIR'bam/derived_'$i'1.fq' $OUT_DIR'fastq/'
  mv $OUT_DIR'bam/derived_'$i'2.fq' $OUT_DIR'fastq/'
done"
write_lines(shell_script, "../shell_scripts/art_neutral_sim.sh")
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_neutral_sim.sh $k /workdir/lcwgs-simulation/neutral_sim_with_replacement/ > '/workdir/lcwgs-simulation/nohups/art_neutral_sim_with_replacement_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
DIR=${2:-/workdir/lcwgs-simulation/neutral_sim/}

OUT_DIR=$DIR'rep_'$REP_ID'/'
N_CORE_MAX=40
## merge
COUNT=0
for k in {1..1000}; do
  i=$((2*$k-1))
  j=$((2*$k))
  ## merge
  samtools merge $OUT_DIR'bam/sample_'$k'.bam' \\
  $OUT_DIR'bam/derived_'$i'.bam' \\
  $OUT_DIR'bam/derived_'$j'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## sort
COUNT=0
for k in {1..1000}; do
  samtools sort -o $OUT_DIR'bam/sample_'$k'_sorted.bam' $OUT_DIR'bam/sample_'$k'.bam' &
  COUNT=$(( COUNT + 1 ))
  if [ $COUNT == $N_CORE_MAX ]; then
    wait
    COUNT=0
  fi
done
wait
## subsample
COUNT=0
for k in {1..1000}; do
  for j in {1,2,4}; do
    samtools view \\
    -s `awk -v j=$j 'BEGIN { print j / 20 }'` \\
    -b $OUT_DIR'bam/sample_'$k'_sorted.bam' \\
    > $OUT_DIR'bam/sample_'$k'_sorted_'$j'x.bam' &
    COUNT=$(( COUNT + 1 ))
    if [ $COUNT == $N_CORE_MAX ]; then
      wait
      COUNT=0
    fi
  done
done
wait
for k in {1..1000}; do
  i=$((2*$k-1))
  j=$((2*$k))
  ## delete intermediate files
  #rm $OUT_DIR'bam/derived_'$i'.bam' 
  #rm $OUT_DIR'bam/derived_'$j'.bam' 
  rm $OUT_DIR'bam/sample_'$k'.bam' 
done"
write_lines(shell_script, "../shell_scripts/merge_sort_subsample_neutral_sim.sh")
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_neutral_sim.sh $k /workdir/lcwgs-simulation/neutral_sim_with_replacement/ > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_neutral_sim_with_replacement_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40,80,160)){
    for (i in 1:sample_size){
      if (i==1){
        write_lines(paste0("/workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_1/bam/sample_", i, "_sorted_", coverage, "x.bam"), paste0("../neutral_sim_with_replacement/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"))
      } else {
        write_lines(paste0("/workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_1/bam/sample_", i, "_sorted_", coverage, "x.bam"), paste0("../neutral_sim_with_replacement/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"), append = T)
      }
    }
  }
}
```