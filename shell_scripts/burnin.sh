#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/sim/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/sim/rep_$REP_ID
fi
# Run SLiM 
/programs/SLiM-3.3/bin/slim 
  -d REP_ID=$REP_ID  
  -d MUTATION_RATE=2e-7 
  -d REC_RATE=1e-8 
  -d CHR_LENGTH=3000 
  -d POP_SIZE=1000 
  -d SAMPLE_SIZE=200 
  -d "OUT_PATH='/workdir/lcwgs-simulation/sim/'" 
  -d "BurninFilename='Burnin.txt'" 
  /workdir/lcwgs-simulation/slim_scripts/burnin.slim
