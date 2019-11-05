#!/bin/bash
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
/programs/SLiM-3.3/bin/slim \
  -d REP_ID=$REP_ID  \
  -d MUTATION_RATE=2e-7 \
  -d REC_RATE=1e-8 \
  -d CHR_LENGTH=30000000 \
  -d POP_SIZE=1000 \
  -d SAMPLE_SIZE=2000 \
  -d "OUT_PATH='/workdir/lcwgs-simulation/neutral_sim/'" \
  /workdir/lcwgs-simulation/slim_scripts/neutral_sim.slim
