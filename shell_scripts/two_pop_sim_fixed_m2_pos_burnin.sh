#!/bin/bash
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
/workdir/programs/SLiM_build/slim \
  -d REP_ID=$REP_ID  \
  -d MUTATION_RATE=100e-8 \
  -d REC_RATE=250e-8 \
  -d MIGRATION_RATE=0.01 \
  -d CHR_LENGTH=30000000 \
  -d POP_SIZE=500 \
  -d "OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'" \
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_burnin.slim
