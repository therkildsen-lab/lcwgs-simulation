#!/bin/bash

## This script is used to build the bow tie reference index. 
# Run this only when working with a new reference that has not been formatted for bowtie2

REFERENCE=$1 # path to reference fasta file and file name, e.g /workdir/cod/reference_seqs/gadMor2.fasta
REFNAME=$2 # reference name to add to output files, e.g. gadMor2

## Build the reference index
REFBASENAME="${REFERENCE%.*}"
bowtie2-build $REFERENCE $REFBASENAME
