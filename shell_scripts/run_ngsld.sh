#!/bin/bash

for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    N_SITE=`zcat /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd/bam_list_${SAMPLE_SIZE}_${COVERAGE}x.pos.gz | wc -l`
    /workdir/programs/ngsLD/ngsLD \
    --geno /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd/bam_list_${SAMPLE_SIZE}_${COVERAGE}x.beagle.gz \
    --pos /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd/bam_list_${SAMPLE_SIZE}_${COVERAGE}x.pos.gz \
    --n_ind $SAMPLE_SIZE \
    --n_sites $((N_SITE - 1)) \
    --out /workdir/lcwgs-simulation/neutral_sim/rep_1/angsd/bam_list_${SAMPLE_SIZE}_${COVERAGE}x.ld \
    --probs \
    --rnd_sample 1 \
    --seed 42 \
    --max_kb_dist 5 \
    --min_maf 0.1 \
    --n_threads 20
    done
done
