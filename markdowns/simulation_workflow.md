Simulation workflow
================

``` r
library(tidyverse)
```

Create a shell script to run SLiM with nohup
--------------------------------------------

``` r
shell_script <- "#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/sim/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/sim/rep_$REP_ID
fi
# Run SLiM 
/programs/SLiM-3.3/bin/slim \\
  -d REP_ID=$REP_ID  \\
  -d MUTATION_RATE=2e-7 \\
  -d REC_RATE=1e-8 \\
  -d CHR_LENGTH=30000000 \\
  -d POP_SIZE=1000 \\
  -d SAMPLE_SIZE=200 \\
  -d \"OUT_PATH='/workdir/lcwgs-simulation/sim/'\" \\
  -d \"BurninFilename='Burnin.txt'\" \\
  /workdir/lcwgs-simulation/slim_scripts/burnin.slim"
write_lines(shell_script, "../shell_scripts/burnin.sh")
```

Run the shell script for SLiM simulation on server
--------------------------------------------------

``` bash
for k in {1..10}; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/burnin.sh $k > '/workdir/lcwgs-simulation/nohups/burnin_'$k'.nohup' &
done
```

Create a shell script to run ART with nohup
-------------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/sim/rep_'$REP_ID'/'
for i in {1..200}; do
  echo '>rep_'$REP_ID > $OUT_DIR'temp_for_art.fasta'
  tail -n 1 'derived_'$i'.fasta' >> $OUT_DIR'temp_for_art.fasta'
  /workdir/programs/art_bin_MountRainier/art_illumina \\
  -ss HS25 \\
  -sam \\
  -i $OUT_DIR'temp_for_art.fasta' \\
  -p \\
  -na \\
  -l 150 \\
  -f 10 \\
  -m 500 \\
  -s 75 \\
  -o 'derived_'$i
samtools view -bS -F 4 $OUT_DIR'derived_'$i'.sam' > $OUT_DIR'derived_'$i'.bam'
rm $OUT_DIR'derived_'$i'.sam'
done"
write_lines(shell_script, "../shell_scripts/run_art.sh")
```

Run the shell script for ART simulation on server
-------------------------------------------------

``` bash
for k in {1..10}; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/run_art.sh $k > '/workdir/lcwgs-simulation/nohups/run_art_'$k'.nohup' &
done
```

Merge, sort, and subsample bam files
------------------------------------

``` r
shell_script <-"#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/sim/rep_'$REP_ID'/'
for k in {1..100}; do
  l=$(($k+100))
  ## merge
  samtools merge $OUTDIR'sample_'$k'.bam' \\
  $OUTDIR'derived_'$k'.bam' \\
  $OUTDIR'derived_'$l'.bam'
  ## sort
  samtools sort -o $OUTDIR'sample_'$k'_sorted.bam' $OUTDIR'sample_'$k'.bam'
  ## subsample
  samtools view \\
  -s 0.05 \\
  -b $OUTDIR'sample_'$k'_sorted.bam' \\
  > $OUTDIR'sample_'$k'_sorted_1x.bam'
  ## delete intermediate files
  #rm $OUTDIR'derived_'$k'.bam' 
  #rm $OUTDIR'derived_'$j'.bam' 
  rm $OUTDIR'sample_'$k'.bam' 
done"
write_lines(shell_script, "../shell_scripts/merge_sort_subsample.sh")
```

Run the shell script for merging, sorting, and subsampling
----------------------------------------------------------

``` bash
for k in {1..10}; do
  nohup bash /workdir/lcwgs-simulation/shell_scripts/merge_sort_subsample.sh $k > '/workdir/lcwgs-simulation/nohups/merge_sort_subsample_'$k'.nohup' &
done
```

Make bam lists
--------------

``` r
for (i in 1:50){
  if (i==1){
    write_lines(paste0("/workdir/lcwgs-simulation/sim/rep_1/sample_", i, "_sorted_1x.bam"), "../sim/rep_1/bam_list_50_1x.txt")
  } else {
    write_lines(paste0("/workdir/lcwgs-simulation/sim/rep_1/sample_", i, "_sorted_1x.bam"), "../sim/rep_1/bam_list_50_1x.txt", append = T)
  }
}
```
