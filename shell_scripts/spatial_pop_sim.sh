#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/spatial_pop_sim/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/spatial_pop_sim/rep_$REP_ID
  cd /workdir/lcwgs-simulation/spatial_pop_sim/rep_$REP_ID
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
  -d MUTATION_RATE=5e-8 \
  -d REC_RATE=1e-8 \
  -d CHR_LENGTH=30000000 \
  -d META_POP_SIDE=3 \
  -d POP_SIZE=500 \
  -d MIGRATION_RATE=0.002 \
  -d SAMPLE_SIZE=160 \
  -d "OUT_PATH='/workdir/lcwgs-simulation/spatial_pop_sim/'" \
  /workdir/lcwgs-simulation/slim_scripts/spatial_pop_sim.slim
