#!/bin/bash
REP_ID=$1
# Create output directory
if [ ! -d /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID ]; then
  mkdir /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID
  cd /workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_$REP_ID
  mkdir angsd
  mkdir bam
  mkdir fasta
  mkdir fastq
  mkdir sample_lists
  mkdir slim
fi
# Copy fasta files
for I in {1..2000}; do
  INDEX=`shuf -i 1-2000 -n 1`
  echo $I'_'$INDEX
  cp '/workdir/lcwgs-simulation/neutral_sim/rep_'$REP_ID'/fasta/derived_'$INDEX'.fasta' '/workdir/lcwgs-simulation/neutral_sim_with_replacement/rep_'$REP_ID'/fasta/derived_'$I'.fasta'
done
