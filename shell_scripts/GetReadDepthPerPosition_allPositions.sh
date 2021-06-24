#!/bin/bash
BAMLIST=$1 # Path to directory of bamfiles to extract depths from (this may be the subsetted bamfiles)
REF=$2 # reference basename to create or find bedfile, e.g. /workdir/Cod/ReferenceSeqs/gadMor2
MEANSTDFILE=$3 # summary file where mean and stdev of each sample is to be written

######Create bedfile from reference (or subsetted reference; this only need to be done once, then comment out)
#samtools faidx $REF'.fasta' #index the fasta file if not done already
#awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $REF'.fasta.fai' > $REF'.bed' # create a bedfile, subset reference bedfile using grep

printf 'Mean\tStDev\n' > $MEANSTDFILE #print header for summary output file
for SAMPLEFILE in `cat $BAMLIST`; do #for each file in the bamlist
####Uncomment next line to save outfiles for each sample
#OUTFILE=`echo "$SAMPLEFILE" | sed -r 's/.bam/.depth/g'` # create an outfile name: depth at each position in the genome will be printed to this file

####Otherwise only save a few files (1 in every ~30)
if test `shuf -i 1-30 -n 1` = 1; then
OUTFILE=`echo "$SAMPLEFILE" | sed -r 's/.bam/.depth/g'` # create an outfile name
else
OUTFILE=Temp.depth #Comment out to generate depth files for each sample
fi

samtools depth -a $SAMPLEFILE | cut -f3 > $OUTFILE #use samtools depth function to get depth, cut out site information to save space
awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {
          printf "%f %f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' $OUTFILE >> $MEANSTDFILE
done

tail -n +2 $MEANSTDFILE | awk '{sumsq+=$2^2; mu+=$1} END {printf "%f %f \n", mu, sqrt(sumsq)}' >> $MEANSTDFILE
rm Temp.depth  #Comment out if saving depth files for each sample

#nohup /workdir/Therkildsen/LowCoverageWorkshop/GetReadDepthPerPosition.sh /workdir/Cod/GoM/SampleLists/GoM_Bamlist.txt /workdir/Cod/ReferenceSeqs/gadMor2_LG12 GoM_ReadDepthPerPosition_LG12.txt >& SamtoolsDepth.nohup &
#nohup /workdir/Therkildsen/LowCoverageWorkshop/GetReadDepthPerPosition.sh /workdir/Cod/Greenland/SampleLists/Greenland_Bamlist_final.txt /workdir/Cod/ReferenceSeqs/gadMor2 Greenland_ReadDepthPerPosition.txt >& SamtoolsDepth.nohup &

