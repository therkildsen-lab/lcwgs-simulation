#!/bin/bash
BASE_DIR='/workdir/lcwgs-simulation/two_pop_sim/rep_1/'
for SAMPLE_SIZE in {5,10,20,40,80,160}; do
  for COVERAGE in {0.25,0.5,1,2,4,8}; do
    # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)				
		/workdir/programs/angsd0.931/angsd/misc/realSFS \
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS'
		# Get alpha_beta
		/workdir/programs/angsd0.931/angsd/misc/realSFS fst index \
    $BASE_DIR'angsd/bam_list_p1_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    $BASE_DIR'angsd/bam_list_p2_'$SAMPLE_SIZE'_'$COVERAGE'x.saf.idx' \
    -sfs $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.2dSFS' \
    -fstout $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta'
		# Get alpha_beta.txt
		/workdir/programs/angsd0.931/angsd/misc/realSFS fst print \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.fst.idx' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt'
		# Get fst
		awk '{ print $0 "	" $3 / $4 }' \
    $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.alpha_beta.txt' \
    > $BASE_DIR'angsd/bam_list_'$SAMPLE_SIZE'_'$COVERAGE'x.fst'
  done
done
