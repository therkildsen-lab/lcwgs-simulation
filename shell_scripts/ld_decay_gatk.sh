#!/bin/bash

for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    echo /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd_gatk/bam_list_${SAMPLE_SIZE}_${COVERAGE}x.ld \
    >> /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd_gatk/bam_list_${SAMPLE_SIZE}.ld.list
  done
  Rscript --vanilla --slave /workdir/programs/ngsLD/scripts/fit_LDdecay.R \
  --ld_files /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd_gatk/bam_list_${SAMPLE_SIZE}.ld.list \
  --out /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd_gatk/ld_decay_${SAMPLE_SIZE}.pdf \
  --fit_level 2 \
  --n_ind $SAMPLE_SIZE \
  --fit_boot 1000 \
  --col 5 \
  > /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd_gatk/ld_decay_${SAMPLE_SIZE}.out &
done

    
