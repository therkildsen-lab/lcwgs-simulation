Simulation workflow for neutral simulation
================

-   [Neutral simulation](#neutral-simulation)
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
    -   [Get shell script for estimating thetas and Tajima's D](#get-shell-script-for-estimating-thetas-and-tajimas-d)

``` r
library(tidyverse)
```

Neutral simulation
==================

Create a shell script to run SLiM with nohup
--------------------------------------------

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/neutral_sim/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/neutral_sim/rep_$REP_ID
  cd /workdir/lcwgs-simulation/neutral_sim/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Run SLiM 
/programs/SLiM-3.3/bin/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=2e-7 \\
  -d REC_RATE=1e-8 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=2000 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/neutral_sim/'\" \\
  /workdir/lcwgs-simulation/slim_scripts/neutral_sim.slim"
write_lines(shell_script, "../shell_scripts/neutral_sim.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in 1; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/neutral_sim.sh $k > '/workdir/lcwgs-simulation/nohups/neutral_sim_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
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
  nohup bash /workdir/lcwgs-simulation/shell_scripts/art_neutral_sim.sh $k > '/workdir/lcwgs-simulation/nohups/art_neutral_sim_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
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
  for j in {0.25,0.5,1,2,4,8}; do
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
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample_neutral_sim.sh $k > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_neutral_sim_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40,80,160)){
    for (i in 1:sample_size){
      if (i==1){
        write_lines(paste0("/workdir/lcwgs-simulation/neutral_sim/rep_1/bam/sample_", i, "_sorted_", coverage, "x.bam"), paste0("../neutral_sim/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"))
      } else {
        write_lines(paste0("/workdir/lcwgs-simulation/neutral_sim/rep_1/bam/sample_", i, "_sorted_", coverage, "x.bam"), paste0("../neutral_sim/rep_1/sample_lists/bam_list_", sample_size, "_", coverage, "x.txt"), append = T)
      }
    }
  }
}
```

Index the ancestral fasta file
------------------------------

``` bash
for REP_ID in 1; do
  OUT_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/'
  samtools faidx $OUT_DIR'slim/ancestral.fasta'
done
```

Get shell script for SNP calling
--------------------------------

``` r
## Note that I used -doMajorMinor 5, using the ancestral sequence to determine major and minor alleles
shell_script <-"#!/bin/bash
SAMPLE_SIZE=$1
COVERAGE=$2
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
## SNP calling
/workdir/programs/angsd0.931/angsd/angsd \\
-b $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \\
-anc $BASE_DIR'slim/ancestral.fasta' \\
-out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \\
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \\
-doCounts 1 -doDepth 1 -dumpCounts 1 \\
-P 4 -SNP_pval 1e-6 -rmTriallelic 1e-6 \\
-setMinDepth 2 -minInd 1 -minMaf 0.0005 -minQ 20 \\
>& '/workdir/lcwgs-simulation/nohups/snp_calling_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.log'
## estimate SFS
/workdir/programs/angsd0.931/angsd/misc/realSFS $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \\
> $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs'"
write_lines(shell_script, "../shell_scripts/snp_calling_neutral_sim.sh")
```

Note: To deal with triallelic loci, I tried `-doMaf 2`, but it turned out that `-doMaf 2` does not bin all alleles other than the ancestral allele. It takes a weighted average approach instead, and makes the result more difficult to process. So we shoud still use`-doMaf 2`, with a `-rmTriallelic 1e-6` filter to remove triallelic loci.

Run the shell script for SNP calling
------------------------------------

``` bash
for coverage in {0.25,0.5,1,2,4,8}; do
  for sample_size in {5,10,20,40,80,160}; do
    nohup bash /workdir/lcwgs-simulation/shell_scripts/snp_calling_neutral_sim.sh \
    $sample_size $coverage \
    > '/workdir/lcwgs-simulation/nohups/snp_calling_neutral_sim_1_bam_list_'$sample_size'_'$coverage'x.nohup' &
  done
done
```

Get shell script for estimating thetas and Tajima's D
-----------------------------------------------------

``` bash
## Get saf file
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/angsd \
    -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -doSaf 1 \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -GL 1 \
    -P 1 \
    > '/workdir/lcwgs-simulation/nohups/get_saf_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' &
  done
done
## Get SFS from saf (this need to be broken up into two separate runs due to memory limiations)
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/realSFS \
      $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
      -P 4 \
      >& $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
      2> '/workdir/lcwgs-simulation/nohups/get_sfs_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' < /dev/null &
  done
done
## Estimate theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/angsd \
    -bam $BASE_DIR'sample_lists/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
    -out $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x' \
    -doThetas 1 \
    -doSaf 1 \
    -pest $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.sfs' \
    -anc $BASE_DIR'slim/ancestral.fasta' \
    -GL 1 \
    -P 1 \
    > '/workdir/lcwgs-simulation/nohups/estimate_theta_neutral_sim_1_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.nohup' &
  done
done
## Print per-SNP theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80,160}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat print \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.tsv' &
  done
done
## Do fixed window theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.windowed_thetas.log' &
  done
done
## Do per-chromosome average theta
BASE_DIR='/workdir/lcwgs-simulation/neutral_sim/rep_1/'
for COVERAGE in {0.25,0.5,1,2,4,8}; do
  for SAMPLE_SIZE in {5,10,20,40,80}; do
    nohup /workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.thetas.idx' \
    -outnames $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.average_thetas.log' &
  done
done
```
