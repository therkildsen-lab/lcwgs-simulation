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
  -d MUTATION_RATE=100e-9 \
  -d REC_RATE=250e-9 \
  -d MIGRATION_RATE=0.001 \
  -d CHR_LENGTH=30000000 \
  -d POP_SIZE=5000 \
  -d SAMPLE_SIZE=1000 \
  -d SELECTION_COEFF=0.08 \
  -d N_M2=11 \
  -d M2_FREQUENCY=1 \
  -d "ANCESTRAL_FASTA='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/ancestral.fasta'" \
  -d "BURNIN_FILE='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/burnin_full_output.txt'" \
  -d "OUT_PATH='/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/'" \
  /workdir/lcwgs-simulation/slim_scripts/two_pop_sim_selection.slim
