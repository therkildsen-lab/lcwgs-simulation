#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
N_JOB_MAX=4
COUNT=0
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
      for POP in {1,2}; do
      ## Get MAF
      /workdir/programs/angsd0.931/angsd/angsd \
      -b $BASE_DIR'sample_lists/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -anc $BASE_DIR'slim/ancestral.fasta' \
      -out $BASE_DIR'angsd/bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x' \
      -sites $BASE_DIR'angsd/global_snp_list_bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.txt' \
      -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 5 \
      -doCounts 1 -doDepth 1 -dumpCounts 1 \
      -P 8 \
      -setMinDepth 1 -minInd 1 -minQ 20 \
      >& '/workdir/lcwgs-simulation/nohups/snp_calling_two_pop_sim_bam_list_p'$POP'_'$SAMPLE_SIZE'_'$COVERAGE'x.log' & 
      ## Submit two jobs at a time
      COUNT=$(( COUNT + 1 ))
      if [ $COUNT == $N_JOB_MAX ]; then
       wait
       COUNT=0
      fi
    done
  done
done
