#!/bin/bash
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
/workdir/programs/SLiM_build/slim \
  -d REP_ID=$REP_ID  \
  -d MUTATION_RATE=2e-7 \
  -d SELECTION_COEFF=0.10 \
  -d MIGRATION_RATE=0.05 \
  -d REC_RATE=2e-8 \
  -d CHR_LENGTH=30000000 \
  -d POP_SIZE=1000 \
  -d SAMPLE_SIZE=1000 \
  -d N_M2=6 \
  -d "OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/'" \
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_fixed_m2_pos.slim
