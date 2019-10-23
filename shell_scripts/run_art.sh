#!/bin/bash
REP_ID=$1
OUT_DIR='/workdir/lcwgs-simulation/sim/rep_'$REP_ID'/'
for i in {1..200}; do
  echo '>rep_'$REP_ID > $OUT_DIR'temp_for_art.fasta'
  tail -n 1 'derived_'$i'.fasta' >> $OUT_DIR'temp_for_art.fasta'
  /workdir/programs/art_bin_MountRainier/art_illumina \
  -ss HS25 \
  -sam \
  -i $OUT_DIR'temp_for_art.fasta' \
  -p \
  -na \
  -l 150 \
  -f 10 \
  -m 500 \
  -s 75 \
  -o 'derived_'$i
samtools view -bS -F 4 $OUT_DIR'derived_'$i'.sam' > $OUT_DIR'derived_'$i'.bam'
rm $OUT_DIR'derived_'$i'.sam'
done
