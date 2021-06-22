Assessing the effect of coverage on estimates of Fst and population
structure for empirical data
================
Arne Jacobs

Load R-packages

``` r
library(tidyverse)
library(cowplot)
library(patchwork)
```

-----

# Download whole genome re-sequencing data for Heliconius erato subspecies from NCBI:

Data from: Van Belleghem, S. M., Rastas, P., Papanicolaou, A., Martin,
S. H., Arias, C. F., Supple, M. A., … Papa, R. (2017). Complex modular
architecture around a simple toolkit of wing pattern genes. Nature
Ecology & Evolution, 1(3), 52.

``` bash
nohup fastq-dump --split-files --gzip  SRS1618041 > SRS1618041_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618053 > SRS1618053_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618075 > SRS1618075_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618086 > SRS1618086_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618096 > SRS1618096_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1617994 > SRS1617994_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618005 > SRS1618005_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618012 > SRS1618012_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618013 > SRS1618013_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618014 > SRS1618014_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618015 > SRS1618015_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618016 > SRS1618016_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618017 > SRS1618017_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618018 > SRS1618018_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618032 > SRS1618032_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618033 > SRS1618033_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618034 > SRS1618034_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618044 > SRS1618044_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618045 > SRS1618045_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618046 > SRS1618046_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618047 > SRS1618047_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618057 > SRS1618057_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618056 > SRS1618056_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618058 > SRS1618058_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618059 > SRS1618059_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618060 > SRS1618060_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618061 > SRS1618061_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618062 > SRS1618062_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618063 > SRS1618063_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618065 > SRS1618065_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618066 > SRS1618066_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618067 > SRS1618067_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618068 > SRS1618068_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618069 > SRS1618069_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618070 > SRS1618070_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618071 > SRS1618071_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618072 > SRS1618072_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618073 > SRS1618073_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618074 > SRS1618074_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618076 > SRS1618076_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618077 > SRS1618077_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618082 > SRS1618082_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618083 > SRS1618083_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618084 > SRS1618084_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618085 > SRS1618085_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618087 > SRS1618087_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618088 > SRS1618088_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618089 > SRS1618089_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618090 > SRS1618090_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618091 > SRS1618091_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618092 > SRS1618092_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618093 > SRS1618093_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618094 > SRS1618094_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618095 > SRS1618095_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618098 > SRS1618098_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618100 > SRS1618100_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618101 > SRS1618101_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618102 > SRS1618102_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618103 > SRS1618103_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618107 > SRS1618107_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1617995 > SRS1617995_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1617996 > SRS1617996_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1617997 > SRS1617997_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1617999 > SRS1617999_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618000 > SRS1618000_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618001 > SRS1618001_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618002 > SRS1618002_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618008 > SRS1618008_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618009 > SRS1618009_fastq_dump_nohup.log &
nohup fastq-dump --split-files --gzip  SRS1618010 > SRS1618010_fastq_dump_nohup.log &
```

-----

# Map raw reads to Heliconius erato demophoon reference genome using the low-coverage pipeline::

Build bowtie2 index:

``` bash
nohup sh /workdir/data-processing/scripts/build_bowtie_ref_index.sh /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa.gz helico_ed_v1_ref > helico_ed_v1_ref_buildbowtie_index.sh &
```

Map and process reads using bowtie2:

``` bash
nohup sh /workdir/arne/lcwgs_empirical/heliconius_data/low_coverage_mapping_heli.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heliconius_belleghem2017_samplelist.txt /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heliconius_belleghem2017_sampletable.txt /workdir/arne/lcwgs_empirical/heliconius_data/fastq_raw/ /workdir/arne/lcwgs_empirical/heliconius_data/ _1.fastq.gz _2.fastq.gz very-sensitive /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa.gz helico_ed_v1_ref > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/heliconius_low_cov_mapping.log &
```

Create list of all bam files:

``` bash
ls *sorted.bam  > helico_lcwgs_bamlist.txt
```

Deduplicate and clip overlapping read pairs:

``` bash
nohup sh /workdir/data-processing/scripts/deduplicate_clipoverlap.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist.txt /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heliconius_belleghem2017_bam_sampletable.txt /workdir/arne/lcwgs_empirical/heliconius_data/ helico_ed_v1_ref > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/heliconius_dedup_clip.log &
```

Indel realignment using picard tools:

``` bash
#Create a sequence dictionary
java -Xmx60g -jar /programs/picard-tools-2.9.0/picard.jar CreateSequenceDictionary R=/workdir/arne/pmajor_bcn_wgs/reference_seq/Heliconius_erato_demophoon_v1_-_scaffolds.fa O=/workdir/arne/pmajor_bcn_wgs/reference_seq/Heliconius_erato_demophoon_v1_-_scaffolds.dict

#Perform indel realignment
nohup sh realign_indels.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist.list /workdir/arne/lcwgs_empirical/heliconius_data/ /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa helico_ed_v1_ref > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/heliconius_realign.log &
```

-----

### Estimating depth of coverage and standard deviation for each sample:

``` bash
nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli.nohup &
```

-----

# Downsampling of NGS data, and estimation of genetic differentiation and population structure:

Downsampling of bam files to varying coverages (8x, 4x, 2x, 1x, 0.5x,
0.25x):

``` bash
nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 30000000 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 8x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &

nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 15000000 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 4x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &

nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 7500000 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 2x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &

nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 3250000 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 1x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &

nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 1625000 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 0.5x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &

nohup ./subsample_bam.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamlist_subsampling.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ 812500 /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/helico_lcwgs_bamtable.txt 0.25x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/nohup_helico_lcwgs_bam_subsampling_8x.log &
```

Estimate depth of coverage for each downsampled dataset

``` bash
nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_8x.nohup &

nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_4x.nohup &

nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_2x.nohup &

nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_1x.nohup &

nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_05x.nohup &

nohup sh GetReadDepthPerPosition_allPositions.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heli_readDepthPerPosition_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/bam/ > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/samtools_depth_all_heli_025x.nohup &
```

-----

## Global SNP calling for subsampled datasets on scaffold 1801 (optix region)

### SNP calling for the dataset at 8x coverage

SNP calling

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 1272 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_8x_cov 5 1272 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_8x_logfile.nohup &
```

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_8x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_8x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_8xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_8x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_8x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_8xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_8x_optix_.saf.idx heli_nonbar_8x_optix_.saf.idx > heli_redbar_nonbar_8x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_8x_optix_.saf.idx heli_nonbar_8x_optix_.saf.idx -sfs heli_redbar_nonbar_8x.2DSFS -fstout heli_redbar_nonbar_8xcov > heli_redbar_nonbar_8xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_8xcov.fst.idx > heli_redbar_nonbar_8xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_8xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_8xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_8x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_8xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_8x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_8x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point(shape = 21, colour = "black", fill = "gray90", size=2) +
  geom_smooth(method = "loess", se = F, colour="steelblue", span=0.2) +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
plot.heli_pooled_phen_fst_50kb20kb_8x
```

### 4x coverage

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 636 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_4x_cov 5 636 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_4x_logfile.nohup &
```

Estimate Fst between phenotypes pooled across subspecies:

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_4x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_4x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_4xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_4x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_4x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_4xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_4x_optix_.saf.idx heli_nonbar_4x_optix_.saf.idx > heli_redbar_nonbar_4x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_4x_optix_.saf.idx heli_nonbar_4x_optix_.saf.idx -sfs heli_redbar_nonbar_4x.2DSFS -fstout heli_redbar_nonbar_4xcov > heli_redbar_nonbar_4xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_4xcov.fst.idx > heli_redbar_nonbar_4xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_4xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_4xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_4x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_4xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_4x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_4x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point(shape = 21, colour = "black", fill = "gray90", size=2) +
  geom_smooth(method = "loess", se = F, colour="steelblue", span=0.2) +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
plot.heli_pooled_phen_fst_50kb20kb_4x
```

### 2x coverage

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 318 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_2x_cov 5 318 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_2x_logfile.nohup &

/programs/angsd_20180926/angsd/angsd sites index Global_SNPList_heli_globalsnp_2x_cov.txt
```

-----

#### Estimate Fst between phenotypes pooled across subspecies:

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_2x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_2x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_2xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_2x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_2x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_2xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_2x_optix_.saf.idx heli_nonbar_2x_optix_.saf.idx > heli_redbar_nonbar_2x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_2x_optix_.saf.idx heli_nonbar_2x_optix_.saf.idx -sfs heli_redbar_nonbar_2x.2DSFS -fstout heli_redbar_nonbar_2xcov > heli_redbar_nonbar_2xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_2xcov.fst.idx > heli_redbar_nonbar_2xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_2xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_2xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_2x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_2xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_2x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_2x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point(shape = 21, colour = "black", fill = "gray90", size=2) +
  geom_smooth(method = "loess", se = F, colour="steelblue", span=0.2) +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
plot.heli_pooled_phen_fst_50kb20kb_2x
```

### 1x coverage

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 159 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_1x_cov 5 159 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_1x_logfile.nohup &
```

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_1x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_1x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_1xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_1x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_1x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_1xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_1x_optix_.saf.idx heli_nonbar_1x_optix_.saf.idx > heli_redbar_nonbar_1x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_1x_optix_.saf.idx heli_nonbar_1x_optix_.saf.idx -sfs heli_redbar_nonbar_1x.2DSFS -fstout heli_redbar_nonbar_1xcov > heli_redbar_nonbar_1xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_1xcov.fst.idx > heli_redbar_nonbar_1xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_1xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_1xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_1x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_1xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_1x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_1x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point(shape = 21, colour = "black", fill = "gray90", size=2) +
  geom_smooth(method = "loess", se = F, colour="steelblue", span=0.2) +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
plot.heli_pooled_phen_fst_50kb20kb_1x
```

### 0.5x coverage

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 80 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_05x_cov 5 80 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_05x_logfile.nohup &

/programs/angsd_20180926/angsd/angsd sites index Global_SNPList_heli_globalsnp_05x_cov.txt
```

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_05x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_05x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_05xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_05x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_05x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_05xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_05x_optix_.saf.idx heli_nonbar_05x_optix_.saf.idx > heli_redbar_nonbar_05x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_05x_optix_.saf.idx heli_nonbar_05x_optix_.saf.idx -sfs heli_redbar_nonbar_05x.2DSFS -fstout heli_redbar_nonbar_05xcov > heli_redbar_nonbar_05xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_05xcov.fst.idx > heli_redbar_nonbar_05xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_05xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_05xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_05x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_05xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_05x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_05x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point(shape = 21, colour = "black", fill = "gray90", size=2) +
  geom_smooth(method = "loess", se = F, colour="steelblue", span=0.2) +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
plot.heli_pooled_phen_fst_50kb20kb_05x
```

### 0.25x coverage

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 40 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/ heli_globalsnp_025x_cov 5 40 40 30 Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_025x_logfile.nohup &

/programs/angsd_20180926/angsd/angsd sites index Global_SNPList_heli_globalsnp_025x_cov.txt
```

Generate saf file for each phenotype:

``` bash
nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_redbar_bamlist_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_redbar_025x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_025x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_redbar_025xcov.nohup &

nohup ./angsd_gl_saf_bypop_optix.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/heli_nobar_bamlist_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/genetic_differentiation_optix/ heli_nonbar_025x_optix_ /workdir/arne/lcwgs_empirical/heliconius_data/snp_data/Global_SNPList_heli_globalsnp_025x_cov.txt Herato1801: > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/angsd_pop_saf_heli_nonbar_025xcov.nohup &
```

Generate 2dSFS

``` bash
# Calculate the 2D site-frequency spectrum (extract the spectrum from the output):
nohup /programs/angsd_20180926/angsd/misc/realSFS heli_redbar_025x_optix_.saf.idx heli_redbar_025x_optix_.saf.idx > heli_redbar_nonbar_025x.2DSFS &
```

Estimate Fst

``` bash
# index to analyze same sites:
nohup /programs/angsd_20180926/angsd/misc/realSFS fst index heli_redbar_025x_optix_.saf.idx heli_nonbar_025x_optix_.saf.idx -sfs heli_redbar_nonbar_025x.2DSFS -fstout heli_redbar_nonbar_025xcov > heli_redbar_nonbar_025xcov.index.nohup.log &

#get the global estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats heli_redbar_nonbar_025xcov.fst.idx > heli_redbar_nonbar_025xcov.globalfst &

#get the sliding window estimate
nohup /programs/angsd_20180926/angsd/misc/realSFS fst stats2 heli_redbar_nonbar_025xcov.fst.idx -win 50000 -step 20000 > heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_025xcov.txt &
```

Plot Fst in 50kb sliding windows with 20kb step:

``` r
heli_pooled_phen_fst_50kb20kb_025x <- read_delim("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_redbar_nonbar_50kbwin_20kbstep_fst_slidingwindow_025xcov.txt", delim = "\t", col_names = T)

plot.heli_pooled_phen_fst_50kb20kb_025x <- ggplot(data = heli_pooled_phen_fst_50kb20kb_025x, aes(x = (midPos/1000000), y=Fst)) +
  geom_point() +
  geom_vline(xintercept = 1250588/1000000, linetype = "dashed") +
  geom_vline(xintercept = 1251208/1000000, linetype = "dashed") +
  xlab("Position along scaffold 1801 [in Mb]") +
  theme_cowplot() 
```

### Gridplot for all 50kb Fst manhattan plots (8x - 0.25x)

``` r
Fst_optix_win.grid <- plot_grid(plot.heli_pooled_phen_fst_50kb20kb_8x,
                          plot.heli_pooled_phen_fst_50kb20kb_4x,
                          plot.heli_pooled_phen_fst_50kb20kb_2x,
                          plot.heli_pooled_phen_fst_50kb20kb_1x,
                          plot.heli_pooled_phen_fst_50kb20kb_05x,
                          ncol = 1)
save_plot(plot = Fst_optix_win.grid, filename = "/Users/arnejacobs/Dropbox/GreatTits_Fst_gridplot_200kbwin_50kbstep_test.pdf", base_height = 10, base_width = 6)
```

-----

### Correlation of Fst values

``` r
cor.test(heli_pooled_phen_fst_50kb20kb_8x$Fst, heli_pooled_phen_fst_50kb20kb_4x$Fst, method = "spearman", exact = F)
cor.test(heli_pooled_phen_fst_50kb20kb_8x$Fst, heli_pooled_phen_fst_50kb20kb_2x$Fst, method = "spearman", exact = F)
cor.test(heli_pooled_phen_fst_50kb20kb_8x$Fst, heli_pooled_phen_fst_50kb20kb_1x$Fst, method = "spearman", exact = F)

cor.test(heli_pooled_phen_fst_50kb20kb_4x$Fst, heli_pooled_phen_fst_50kb20kb_2x$Fst, method = "spearman", exact = F)
cor.test(heli_pooled_phen_fst_50kb20kb_4x$Fst, heli_pooled_phen_fst_50kb20kb_1x$Fst, method = "spearman", exact = F)

cor.test(heli_pooled_phen_fst_50kb20kb_2x$Fst, heli_pooled_phen_fst_50kb20kb_1x$Fst, method = "spearman", exact = F)

df1=as_tibble(merge(heli_pooled_phen_fst_50kb20kb_8x, heli_pooled_phen_fst_50kb20kb_05x, by="midPos"))
cor.test(df1$Fst.x, df1$Fst.y, method = "spearman", exact = F)

df1=as_tibble(merge(heli_pooled_phen_fst_50kb20kb_4x, heli_pooled_phen_fst_50kb20kb_05x, by="midPos"))
cor.test(df1$Fst.x, df1$Fst.y, method = "spearman", exact = F)

df1=as_tibble(merge(heli_pooled_phen_fst_50kb20kb_2x, heli_pooled_phen_fst_50kb20kb_05x, by="midPos"))
cor.test(df1$Fst.x, df1$Fst.y, method = "spearman", exact = F)

df1=as_tibble(merge(heli_pooled_phen_fst_50kb20kb_1x, heli_pooled_phen_fst_50kb20kb_05x, by="midPos"))
cor.test(df1$Fst.x, df1$Fst.y, method = "spearman", exact = F)
```

-----

# Genome-wide Principal Components Analysis (PCA) for heliconius data for different coverages:

## Create global and population specific SNP files for population genomic analyses:

### Call SNPs for each subsampled dataset:

#### 8X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 1272 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_8x_cov 5 1272 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_global_snpcalling_8x_logfile.nohup &
```

#### 4X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 636 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_4x_cov 5 636 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_4x_logfile.nohup &
```

#### 2X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 318 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_2x_cov 5 318 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_2x_logfile.nohup &
```

#### 1X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 159 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_1x_cov 5 159 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_1x_logfile.nohup &
```

#### 0.5X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 80 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_05x_cov 5 80 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_05x_logfile.nohup &
```

#### 0.25X

Used filters: \* minCov: 0.1 x N.ind = 5 \* maxCov: cov*N.ind + (2 *
cov*N.ind) = 40 * minInd: 75% = 40 \* minQ = 30

``` bash
nohup sh Angsd_GlobalSNPcalling.sh /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_025x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ heli_globalsnp_025x_cov 5 40 40 30 > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_snpcalling_025x_logfile.nohup &
```

## Principal components analysis for each coverage

Estimation of covariance matrices using ANGSD

``` bash
nohup sh ./angsd_getIBS.sh /workdir/arne/lcwgs_empirical/heliconius_data/bam/ /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_8x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/Global_SNPList_heli_globalsnp_8x_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ all_filtsnps_IBS_8x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_IBS_8x_logfile.nohup &

nohup sh ./angsd_getIBS.sh /workdir/arne/lcwgs_empirical/heliconius_data/bam/ /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_4x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/Global_SNPList_heli_globalsnp_4x_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ all_filtsnps_IBS_4x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_IBS_4x_logfile.nohup &

nohup sh ./angsd_getIBS.sh /workdir/arne/lcwgs_empirical/heliconius_data/bam/ /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_2x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/Global_SNPList_heli_globalsnp_2x_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ all_filtsnps_IBS_2x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_IBS_2x_logfile.nohup &

nohup sh ./angsd_getIBS.sh /workdir/arne/lcwgs_empirical/heliconius_data/bam/ /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_1x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/Global_SNPList_heli_globalsnp_1x_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ all_filtsnps_IBS_1x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_IBS_1x_logfile.nohup &

nohup sh ./angsd_getIBS.sh /workdir/arne/lcwgs_empirical/heliconius_data/bam/ /workdir/arne/lcwgs_empirical/heliconius_data/sampleinfo/bamlist_05x.txt /workdir/arne/lcwgs_empirical/heliconius_data/Heliconius_erato_demophoon_v1_-_scaffolds.fa /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/Global_SNPList_heli_globalsnp_05x_cov.txt /workdir/arne/lcwgs_empirical/heliconius_data/population_structure/ all_filtsnps_IBS_05x > /workdir/arne/lcwgs_empirical/heliconius_data/logfiles/output_heli_IBS_05x_logfile.nohup &
```

### PCA and plotting

8x coverage

``` r
cov.8x.m=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/all_filtsnps_IBS_8x.covMat"))
ind.8x.df=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/bamlist_8x.txt"))
cov8x.pca=eigen(cov.8x.m)
cov8x.pca=as_tibble(cov8x.pca$vectors)
names(cov8x.pca)[1]=c("PC1")
names(cov8x.pca)[2]=c("PC2")
cov8x.pca.df=cbind(ind.8x.df,cov8x.pca)
PCA8x=ggplot(data=cov8x.pca.df, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill = V1), size = 4, shape=21) +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  scale_fill_viridis_d() +
  theme_cowplot()+
  theme(legend.position = "NA")
PCA8x
```

4x coverage

``` r
cov.4x.m=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/all_filtsnps_IBS_4x.covMat"))
ind.4x.df=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/bamlist_4x.txt"))
cov4x.pca=eigen(cov.4x.m)
cov4x.pca=as_tibble(cov4x.pca$vectors)
names(cov4x.pca)[1]=c("PC1")
names(cov4x.pca)[2]=c("PC2")
cov4x.pca.df=cbind(ind.4x.df,cov4x.pca)
PCA4x=ggplot(data=cov4x.pca.df, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill = V1), size = 4, shape=21) +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  scale_fill_viridis_d() +
  theme_cowplot()+
  theme(legend.position = "NA")
PCA4x
```

2x coverage

``` r
cov.2x.m=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/all_filtsnps_IBS_2x.covMat"))
ind.2x.df=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/bamlist_2x.txt"))
cov2x.pca=eigen(cov.2x.m)
cov2x.pca=as_tibble(cov2x.pca$vectors)
names(cov2x.pca)[1]=c("PC1")
names(cov2x.pca)[2]=c("PC2")
cov2x.pca.df=cbind(ind.2x.df,cov2x.pca)
PCA2x=ggplot(data=cov2x.pca.df, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill = V1), size = 4, shape=21) +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  scale_fill_viridis_d() +
  theme_cowplot()+
  theme(legend.position = "NA")
PCA2x
```

1x coverage

``` r
cov.1x.m=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/all_filtsnps_IBS_1x.covMat"))
ind.1x.df=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/bamlist_1x.txt"))
cov1x.pca=eigen(cov.1x.m)
cov1x.pca=as_tibble(cov1x.pca$vectors)
names(cov1x.pca)[1]=c("PC1")
names(cov1x.pca)[2]=c("PC2")
cov1x.pca.df=cbind(ind.1x.df,cov1x.pca)
PCA1x=ggplot(data=cov1x.pca.df, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill = V1), size = 4, shape=21) +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  scale_fill_viridis_d() +
  theme_cowplot()+
  theme(legend.position = "NA")
PCA1x
```

0.5x coverage

``` r
cov.05x.m=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/all_filtsnps_IBS_05x.covMat"))
ind.05x.df=as.matrix(read.table("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/bamlist_05x.txt"))
cov05x.pca=eigen(cov.05x.m)
cov05x.pca=as_tibble(cov05x.pca$vectors)
names(cov05x.pca)[1]=c("PC1")
names(cov05x.pca)[2]=c("PC2")
cov05x.pca.df=cbind(ind.05x.df,cov05x.pca)
PCA05x=ggplot(data=cov05x.pca.df, aes(x=PC1,y=PC2)) +
  geom_point(aes(fill = V1), size = 4, shape=21) +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  scale_fill_viridis_d() +
  theme_cowplot()+
  theme(legend.position = "NA")
PCA05x
```

Create grid plot:

``` r
pca.grid= plot_grid(PCA8x, PCA4x, PCA2x, PCA1x, PCA05x, nrow = 5, ncol = 1 ,labels = c("8x","4x","2x","1x","0.5x"))
save_plot(plot = pca.grid, filename = "/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/heli_pca_grid.pdf", base_height = 10, base_width = 3)
```

-----

# Coverage vs Number of SNPs

Correlation betwen the coverage and number of called SNPs (Quadratic
model)

``` r
snp.cov.df=read_csv("/Users/arnejacobs/Dropbox/Cornell_Postdoc/analysis/lcwgs_guide_empirical/analysis/HeliLowCov_Cov_vs_Nsnps.csv")

snp.cov.df$NSNPs2 <- (snp.cov.df$NSNPs)^2
snp.cov.df$Coverage2 <- (snp.cov.df$Coverage)^2
summary(lm(NSNPs2 ~ Coverage, data = snp.cov.df))
```

Plot correlation

``` r
quad.plot=ggplot(data = snp.cov.df, aes(x=Coverage, y=NSNPs)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), size = 1, colour = "gray20") +
  geom_point(size = 4, shape = 21, fill="white", stroke=2) +
  annotate(geom="text", x=6, y=0.01, label="R2 = 0.98, P = 0.00099",color="black") +
  ylab("Number of SNPs") +
  xlab("Depth of coverage") +
  theme_cowplot()
```

-----

# Generate multi-panel plot across analyses using patchwork

``` r
plot = (quad.plot/plot.heli_pooled_phen_fst_50kb20kb_8x/plot.heli_pooled_phen_fst_50kb20kb_1x/plot.heli_pooled_phen_fst_50kb20kb_05x) | (PCA8x / PCA1x / PCA05x)
```
