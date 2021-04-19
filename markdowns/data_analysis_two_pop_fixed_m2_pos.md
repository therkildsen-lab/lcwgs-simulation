Data analysis with simulation of divergent selection on two populations
================

  - [Define all relevant functions](#define-all-relevant-functions)
  - [Standard model (Ne\~50,000 in each
    population)](#standard-model-ne50000-in-each-population)
      - [The model](#the-model)
          - [Read in the ancestral
            states](#read-in-the-ancestral-states)
          - [Read mutation and substitution
            file](#read-mutation-and-substitution-file)
          - [Plot Fst](#plot-fst)
          - [Plot Fst from the Greenland cod project as a
            comparison](#plot-fst-from-the-greenland-cod-project-as-a-comparison)
      - [Inference with Samtool’s GL
        model](#inference-with-samtools-gl-model)
          - [Read in read depth and estimated
            Fst](#read-in-read-depth-and-estimated-fst)
          - [Plot genome-wide average Fst (with no minimum individual
            filter)](#plot-genome-wide-average-fst-with-no-minimum-individual-filter)
          - [Plot the estimated per-SNP Fst (with no minimum individual
            filter)](#plot-the-estimated-per-snp-fst-with-no-minimum-individual-filter)
          - [Plot genome-wide average Fst (with minimum individual
            filter)](#plot-genome-wide-average-fst-with-minimum-individual-filter)
          - [Plot the estimated per-SNP Fst (with minimum individual
            filter)](#plot-the-estimated-per-snp-fst-with-minimum-individual-filter)
          - [Compute and plot the estimated windowed Fst (with no
            minimum individual filter and 1,000bp fixed
            windows)](#compute-and-plot-the-estimated-windowed-fst-with-no-minimum-individual-filter-and-1000bp-fixed-windows)
          - [Selection scan using
            PCAngsd](#selection-scan-using-pcangsd)
      - [Inference with GATK’s GL model](#inference-with-gatks-gl-model)
          - [Fst in 1000bp windows](#fst-in-1000bp-windows)
          - [Selection scan using neutrality test
            stats](#selection-scan-using-neutrality-test-stats)
              - [Read in the data](#read-in-the-data)
              - [Genome wide stats](#genome-wide-stats)
              - [Window-based stats](#window-based-stats)
              - [Tajima’s D](#tajimas-d)
              - [Fay and Wu’s H](#fay-and-wus-h)
      - [RAD seq simulaiton and
        inference](#rad-seq-simulaiton-and-inference)
          - [Get true sample allele
            count](#get-true-sample-allele-count)
          - [Get the true SFS and theta estimators with population
            samples in
            p1](#get-the-true-sfs-and-theta-estimators-with-population-samples-in-p1)
          - [Get sample true MAF and Fst from allele
            counts](#get-sample-true-maf-and-fst-from-allele-counts)
          - [Chromosome-wide stats](#chromosome-wide-stats)
              - [Plot SFS in p1 from RAD
                data](#plot-sfs-in-p1-from-rad-data)
              - [Tajima’s estimator](#tajimas-estimator-1)
              - [Watterson’s estimator](#wattersons-estimator-1)
          - [Plot per SNP Fst](#plot-per-snp-fst)
  - [Two populations with divergent selection, with smaller population
    size ( Ne\~10,000 in each
    population)](#two-populations-with-divergent-selection-with-smaller-population-size-ne10000-in-each-population)
      - [The model](#the-model-1)
          - [Read in the ancestral
            states](#read-in-the-ancestral-states-1)
          - [Read mutation and substitution
            file](#read-mutation-and-substitution-file-1)
          - [Plot Fst](#plot-fst-1)
      - [Inference based on Samtool’s GL
        model](#inference-based-on-samtools-gl-model)
          - [Read in read depth and estimated
            Fst](#read-in-read-depth-and-estimated-fst-1)
          - [Plot the estimated per-SNP Fst (with no minimum individual
            filter)](#plot-the-estimated-per-snp-fst-with-no-minimum-individual-filter-1)
          - [Compute and plot the estimated windowed Fst (with no
            minimum individual filter and 5,000bp fixed
            windows)](#compute-and-plot-the-estimated-windowed-fst-with-no-minimum-individual-filter-and-5000bp-fixed-windows)
          - [Selection scan using
            PCAngsd](#selection-scan-using-pcangsd-1)
      - [Inference based on GATK’s GL
        model](#inference-based-on-gatks-gl-model)
          - [Selection scan using neutrality test
            stats](#selection-scan-using-neutrality-test-stats-1)
              - [Read in the data](#read-in-the-data-1)
              - [Genome wide stats](#genome-wide-stats-1)
              - [Window-based stats](#window-based-stats-1)
      - [RAD seq simulation and
        inference](#rad-seq-simulation-and-inference)
          - [Get true sample allele
            count](#get-true-sample-allele-count-1)
          - [Get the true SFS and theta estimators with population
            samples in
            p1](#get-the-true-sfs-and-theta-estimators-with-population-samples-in-p1-1)
          - [Get sample true MAF and Fst from allele
            counts](#get-sample-true-maf-and-fst-from-allele-counts-1)
          - [Plot SFS in p1 from RAD
            data](#plot-sfs-in-p1-from-rad-data-1)
          - [Get theta estimation in p1 from RAD
            data](#get-theta-estimation-in-p1-from-rad-data)
              - [Tajima’s estimator](#tajimas-estimator-3)
              - [Watterson’s estimator](#wattersons-estimator-3)
          - [Plot per SNP Fst](#plot-per-snp-fst-1)

``` r
library(tidyverse)
library(cowplot)
library(knitr)
library(data.table)
library(RcppCNPy)
```

# Define all relevant functions

``` r
get_ancestral <- function(x){
  read_csv(paste0(x,"slim/ancestral.fasta"), col_types = cols())[[1]] %>%
    str_split(pattern="") %>%
    .[[1]] %>%
    bind_cols(ancestral=., position=1:30000000)
}

get_mutations <- function(x){
  ## Read in the mutation file outputted by SLiM
  mutations <- read_delim(paste0(x, "slim/mutations.txt"), delim = " ", col_names = F, col_types = cols()) %>%
    transmute(population=X4, type=X6, position=X7+1, base=X13, frequency=X12/10000) %>%
    left_join(ancestral, by="position") %>%
    group_by(population, type, position, ancestral, base) %>%
    summarise(frequency=sum(frequency)) %>%
    ungroup()
  ## Read in the substitutions file outputted by SLiM
  ## This is necessary because mutations can happen again after one fixation, so frequencies from the mutation file do not always reflect the true derived allele frequency
  substitutions <- read_delim(paste0(x,"slim/substitutions.txt"), delim = " ", skip=2, col_names = F, col_types = cols()) %>%
    transmute(type=X3, position=X4+1, base=X10, generation=X9, p1=1, p2=1) %>%
    group_by(type, position) %>%
    filter(generation==max(generation)) %>%
    ungroup() %>%
    left_join(ancestral, by="position") %>%
    dplyr::select(-generation) %>%
    filter(base!=ancestral) %>%
    gather(key=population, value=frequency, 4:5) %>%
    arrange(position)
  ## The following steps are necessary because there are complications such as back mutations and triallelic loci in the mutation file
  ## Join mutations and substitutions in a temp table
  mutations_final_temp <-  mutations %>%
    spread(key = base, value=frequency) %>%
    full_join(substitutions, by=c("position", "type", "ancestral", "population")) %>%
    arrange(position) %>%
    mutate(base=ifelse(is.na(base), ancestral, base)) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    mutate(frequency=1-`A` -`C` -`G` -`T`)
  ## More wrangling
  mutations_final <- mutations_final_temp[1:8] %>%
    gather(key=base, value=frequency, 5:8) %>%
    bind_rows(mutations_final_temp[c(1:4, 9:10)]) %>%
    mutate(frequency=ifelse(base==ancestral, 0, frequency)) %>%
    group_by(population, type, position, ancestral) %>%
    summarise(frequency=sum(frequency)) %>%
    ungroup() %>%
    spread(key=population, value=frequency) %>%
    mutate_all(~replace(., is.na(.), 0)) %>% 
    filter(!(p1==1 & p2==1), !(p1==0 & p2==0)) %>%
    mutate(frequency_mean = (p1 + p2)/2, h_t=2*frequency_mean*(1-frequency_mean), h_s=p1*(1-p1) + p2*(1-p2), fst=1-h_s/h_t)
  return(mutations_final)
}

get_estimated_fst <- function(x){
  i=1
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80, 160)){
      ## read in estimated fst
      fst <- read_tsv(paste0(x, "angsd/bam_list_", sample_size, "_", coverage, "x.fst"), col_names = F, col_types = cols()) %>%
        transmute(position=X2, alpha=X3, beta=X4, fst=X5, coverage=coverage, sample_size=sample_size)
      ## read per population depth
      p1_n_ind <- read_tsv(paste0(x, "angsd/bam_list_p1_", sample_size, "_", coverage, "x.mafs.gz"), col_types = cols()) %>%
        transmute(position=position, p1_n_ind=nInd)
      p2_n_ind <- read_tsv(paste0(x, "angsd/bam_list_p2_", sample_size, "_", coverage, "x.mafs.gz"), col_types = cols()) %>%
        transmute(position=position, p2_n_ind=nInd)
      ## join fst with depth_files
      fst_n_ind <- left_join(fst, p1_n_ind, by="position") %>%
        left_join(p2_n_ind, by="position")
      ## compile the final files for plotting
      if (i==1){
        fst_n_ind_final <- fst_n_ind
      } else {
        fst_n_ind_final <- bind_rows(fst_n_ind_final, fst_n_ind)
      }
      i=i+1
    }
  }
  return(fst_n_ind_final)
}

get_estimated_windowed_fst_gatk <- function(x, window_length){
  i=1
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80, 160)){
      ## read in estimated fst
      fst <- read_tsv(paste0(x, "angsd_gatk/bam_list_", sample_size, "_", coverage, "x.alpha_beta.txt"), col_names = F, col_types = cols()) %>%
        transmute(position=X2, alpha=X3, beta=X4, coverage=coverage, sample_size=sample_size) %>%
        mutate(position=cut(position, 
                         breaks=seq(0,40*10^6,window_length),
                         labels=seq(window_length/2,40*10^6-window_length/2,window_length))) %>%
        group_by(position, coverage, sample_size) %>%
        summarise(alpha = sum(alpha), beta = sum(beta), fst = alpha / beta) %>%
        ungroup() %>%
        mutate(position=as.numeric(as.character(position)))
      ## compile the final files for plotting
      if (i==1){
        fst_final <- fst
      } else {
        fst_final <- bind_rows(fst_final, fst)
      }
      i=i+1
    }
  }
  return(fst_final)
}

fixed_windowed_fst <- function(x, window_length){
  mutate(x, position=cut(position, 
                         breaks=seq(0,40*10^6,window_length),
                         labels=seq(window_length/2,40*10^6-window_length/2,window_length))) %>%
    group_by(position, coverage, sample_size) %>%
    summarise(fst=sum(alpha)/sum(beta)) %>%
    ungroup() %>%
    mutate(position=as.numeric(as.character(position)))
}

get_selection_scan <- function(x){
  i=1
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80, 160)){
      ## read in estimated fst
      genome_selection <- npyLoad(paste0(x, "angsd/pcagnsd_bam_list_", sample_size, "_", coverage, "x.selection.npy"))
      genome_selection_sites <- read_table(paste0(x, "angsd/pcagnsd_bam_list_", sample_size, "_", coverage, "x.sites"), col_names = F) %>%
        transmute(coverage = coverage, 
                  sample_size = sample_size,
                  pos = parse_integer(str_remove(X1, "rep_1_")), 
                  chi_squared = genome_selection[,1],
                  neg_log_p_value = -log(1-pchisq(chi_squared, df=1)))
      ## compile the final files for plotting
      if (i==1){
        genome_selection_sites_final <- genome_selection_sites
      } else {
        genome_selection_sites_final <- bind_rows(genome_selection_sites_final, genome_selection_sites)
      }
      i=i+1
    }
  }
  return(genome_selection_sites_final)
}
get_neutrality_stats <- function(x){
  i=1
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80, 160)){
      ## read in estimated fst
      neutrality_stats <- read_tsv(paste0(x, "angsd_gatk/bam_list_p1_", sample_size, "_", coverage, "x.windowed_thetas.idx.pestPG")) %>%
        janitor::clean_names() %>%
        dplyr::select(-1, -2) %>%
        mutate(coverage = coverage, 
                  sample_size = sample_size) %>%
        rename(pos = win_center)
      ## compile the final files for plotting
      if (i==1){
        neutrality_stats_final <- neutrality_stats
      } else {
        neutrality_stats_final <- bind_rows(neutrality_stats_final, neutrality_stats)
      }
      i=i+1
    }
  }
  return(neutrality_stats_final)
}
count_to_maf <- function(ancestral_allele, totA, totC, totG, totT){
  if(ancestral_allele == "A"){
    minor_allele_count <- max(totC, totG, totT)
  } else if(ancestral_allele == "C"){
    minor_allele_count <- max(totA, totG, totT)
  } else if(ancestral_allele == "G"){
    minor_allele_count <- max(totA, totC, totT)
  } else if(ancestral_allele == "T"){
    minor_allele_count <- max(totA, totC, totG)
  }
  maf <- minor_allele_count/sum(totA, totC, totG, totT)
  return(maf)
}
get_sample_allele_count_per_pop <- function(x){
  for (p in 1:2){
    i <- 1
    for (sample_id in 1:160){
      for (genome in 1:2){
        sequence <- read_csv(paste0(x, "fasta/p", p, "_derived_", sample_id, "_", genome, ".fasta"), col_types = cols())[[1]] %>%
          str_split(pattern="") %>%
          .[[1]] %>%
          tibble(base=., position=1:30000000, ancestral=ancestral$ancestral) %>%
          filter(position %in% mutations_final$position)
        allele_count <- transmute(sequence,
                                  A_count = ifelse(base=="A", 1, 0),
                                  C_count = ifelse(base=="C", 1, 0),
                                  G_count = ifelse(base=="G", 1, 0),
                                  T_count = ifelse(base=="T", 1, 0))
        if (i==1){
          allele_count_final <- allele_count
        } else {
          allele_count_final <- allele_count + allele_count_final
        }
        i <- i+1
      }
      if (sample_id %in% c(5,10,20,40,80,160)){
        write_tsv(bind_cols(dplyr::select(sequence, -base), allele_count_final), paste0(x,"slim/p", p, "_", sample_id, "_base_count.tsv")) 
      }
    }
  }
}
allele_count_to_fst <- function(x){
  i <- 1
  for (sample_size in c(5,10,20,40,80,160)){
    base_count_p1 <- read_tsv(paste0(x,"slim/p1_", sample_size, "_base_count.tsv"), col_types = cols())
    maf_p1 <- base_count_p1 %>%
      rowwise() %>%
      transmute(maf = count_to_maf(ancestral, A_count, C_count, G_count, T_count), position=position) %>%
      ungroup()
    base_count_p2 <- read_tsv(paste0(x, "slim/p2_", sample_size, "_base_count.tsv"), col_types = cols())
    maf_p2 <- base_count_p2 %>%
      rowwise() %>%
      transmute(maf = count_to_maf(ancestral, A_count, C_count, G_count, T_count), position=position) %>%
      ungroup()
    fst <- tibble(position=maf_p1$position, p1=maf_p1$maf, p2=maf_p2$maf) %>%
      rowwise() %>%
      mutate(maf_mean = (p1+p2)/2) %>%
      mutate(h_t=2*maf_mean*(1-maf_mean), h_s=p1*(1-p1) + p2*(1-p2), fst=1-h_s/h_t, sample_size=sample_size)
    if (i==1) {
      fst_final <- fst
    } else {
      fst_final <- bind_rows(fst_final, fst)
    }
    i <- i+1
  }
  return(fst_final)
}
subset_fst <- function(n_rad_tag, combined_fst){
  set.seed(1)
  rad_intervals <- sample(1:(30000000-150), n_rad_tag) %>%
    tibble(start=., stop=.+150) %>%
    arrange(by=start)
  n_sites <- sum(1:30000000 %inrange% as.list(rad_intervals))
  i <- 1
  for (n_samples in c(5,10,20,40,80,160)){
    fst <- filter(combined_fst, sample_size == n_samples, position %inrange% as.list(rad_intervals)) %>%
      mutate(n_sites = n_sites)
    if (i==1) {
      fst_final <- fst
    } else {
      fst_final <- bind_rows(fst_final, fst)
    }
    i <- i+1
  }
  return(fst_final)
}
```

# Standard model (Ne\~50,000 in each population)

## The model

### Read in the ancestral states

``` r
ancestral <- get_ancestral("../two_pop_sim_fixed_m2_pos/rep_1/")
```

### Read mutation and substitution file

The target theta is \~ 0.004.

``` r
mutations_final <- get_mutations("../two_pop_sim_fixed_m2_pos/rep_1/")
real_theta_t_p1 <- sum(2*mutations_final$p1*(1-mutations_final$p1))/30000000
real_theta_t_p2 <- sum(2*mutations_final$p2*(1-mutations_final$p2))/30000000
real_theta_t_combined <- sum(2*mutations_final$frequency_mean*(1-mutations_final$frequency_mean))/30000000
real_theta_w_p1 <- filter(mutations_final, p1 > 0, p1 < 1) %>%
  nrow() %>% `/`(30000000*sum(1/(1:(5000*2-1))))
real_theta_w_p2 <- filter(mutations_final, p2 > 0, p2 < 1) %>%
  nrow() %>% `/`(30000000*sum(1/(1:(5000*2-1))))
real_theta_w_combined <- nrow(mutations_final) / (30000000*sum(1/(1:(10000*2-1))))
tibble(theta = c("tajima", "watterson"), 
       p1 = c(real_theta_t_p1, real_theta_w_p1), 
       p2 = c(real_theta_t_p2, real_theta_w_p2),
       combined = c(real_theta_t_combined, real_theta_w_combined)) %>% 
  kable()
```

| theta     |        p1 |        p2 |  combined |
| :-------- | --------: | --------: | --------: |
| tajima    | 0.0036054 | 0.0036273 | 0.0036745 |
| watterson | 0.0030680 | 0.0030931 | 0.0038472 |

``` r
## p1 SFS
filter(mutations_final, p1 > 0, p1 < 1) %>%
  mutate(p1 = round(p1, 4)) %>%
  count(p1) %>%
  ggplot(aes(x=p1, y=n)) +
  geom_line() +
  geom_point(size=1) +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
## p2 SFS
filter(mutations_final, p2 > 0, p2 < 1) %>%
  mutate(p2 = round(p2, 4)) %>%
  count(p2) %>%
  ggplot(aes(x=p2, y=n)) +
  geom_line() +
  geom_point(size=1) +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

`theta_w` is much lower than `theta_p` in p1 and p2, presumably because
selection has strongly reduced the effective migration between the two
populations, creating a genome-wide pattern of population contraction
when each population is evaluated separately.

### Plot Fst

The target mean neutral Fst is \~ 0.015 and the target peak Fst is \~
0.6 (based on cod data). The expected neutral Fst is `1/(1+16Nm)` =
0.01234, but because of the selection stage, I expect Fst to be slightly
higher than this theoretical expectation.

I used `mean(h_s)/mean(h_t)` to calculate average Fst. I am not sure if
this estimator is biased though.

``` r
mutations_final_m1 <- filter(mutations_final, type=="m1")
mutations_final_m2 <- filter(mutations_final, type=="m2")
## genome-wide mean fst
summarise(mutations_final, average_fst = 1-mean(h_s)/mean(h_t))
```

    ## # A tibble: 1 x 1
    ##   average_fst
    ##         <dbl>
    ## 1      0.0158

``` r
## "neutral" mean fst (mean fst at both end of the genome that is not strongly affected by linkage with selected regions)
filter(mutations_final, position <= 1000000 | position >= 29000000) %>% summarise(average_fst = 1-mean(h_s)/mean(h_t))
```

    ## # A tibble: 1 x 1
    ##   average_fst
    ##         <dbl>
    ## 1      0.0136

``` r
ggplot(mutations_final_m1, aes(x=position, y=fst, color=type)) +
  geom_point(size=0.02, alpha=0.5) +
  geom_point(data=mutations_final_m2, aes(x=position, y=fst, color=type)) +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
arrange(mutations_final, desc(fst)) %>%
  head()
```

    ## # A tibble: 6 x 9
    ##   type  position ancestral    p1     p2 frequency_mean   h_t    h_s   fst
    ##   <chr>    <dbl> <chr>     <dbl>  <dbl>          <dbl> <dbl>  <dbl> <dbl>
    ## 1 m2    17500001 A         0.918 0.0071          0.462 0.497 0.0827 0.834
    ## 2 m2     7500001 A         0.904 0.0119          0.458 0.496 0.0987 0.801
    ## 3 m1     7504306 T         0.887 0.0117          0.449 0.495 0.112  0.774
    ## 4 m1     7501513 G         0.898 0.0233          0.460 0.497 0.115  0.769
    ## 5 m1    17503032 G         0.862 0.0062          0.434 0.491 0.125  0.746
    ## 6 m1     7502493 C         0.898 0.0453          0.472 0.498 0.134  0.730

``` r
arrange(mutations_final_m2, desc(fst))
```

    ## # A tibble: 8 x 9
    ##   type  position ancestral     p1      p2 frequency_mean    h_t    h_s    fst
    ##   <chr>    <dbl> <chr>      <dbl>   <dbl>          <dbl>  <dbl>  <dbl>  <dbl>
    ## 1 m2    17500001 A         0.918  0.0071          0.462  0.497  0.0827 0.834 
    ## 2 m2     7500001 A         0.904  0.0119          0.458  0.496  0.0987 0.801 
    ## 3 m2     2500001 G         0.841  0.0055          0.423  0.488  0.139  0.715 
    ## 4 m2    22500001 G         0.823  0.0052          0.414  0.485  0.151  0.689 
    ## 5 m2    10000001 G         0.791  0.00930         0.400  0.480  0.174  0.637 
    ## 6 m2    20000001 G         0.762  0.0061          0.384  0.473  0.187  0.604 
    ## 7 m2     5000001 G         0.665  0.0067          0.336  0.446  0.229  0.486 
    ## 8 m2    27500001 C         0.0268 0               0.0134 0.0264 0.0261 0.0136

### Plot Fst from the Greenland cod project as a comparison

``` r
fst <- read_tsv("../../cod/greenland-cod/angsd/popminind2/ILU2011_UUM2010_bam_list_realigned_mindp161_maxdp768_minind97_minq20_popminind2.fst", col_names = F) %>%
  rename(lg=X1, position=X2, alpha=X3, beta=X4, fst=X5)
## LG03
filter(fst, lg=="LG03") %>%
ggplot(aes(x=position, y=fst)) +
  geom_point(size=0.02, alpha=0.5) +
  theme_cowplot()
## LG08
filter(fst, lg=="LG08") %>%
ggplot(aes(x=position, y=fst)) +
  geom_point(size=0.02, alpha=0.5) +
  theme_cowplot()
## LG19
filter(fst, lg=="LG19") %>%
ggplot(aes(x=position, y=fst)) +
  geom_point(size=0.02, alpha=0.5) +
  theme_cowplot()
## Mean Fst at LG19
filter(fst, lg=="LG19") %>% summarise(mean_fst_neutral=sum(alpha)/sum(beta))
```

## Inference with Samtool’s GL model

### Read in read depth and estimated Fst

``` r
fst_n_ind_final <- get_estimated_fst("../two_pop_sim_fixed_m2_pos/rep_1/")
```

### Plot genome-wide average Fst (with no minimum individual filter)

``` r
group_by(fst_n_ind_final, sample_size, coverage) %>%
  count() %>%
  pivot_wider(names_from = sample_size, values_from = n)
average_fst_plot <- fst_n_ind_final %>%
  group_by(coverage, sample_size) %>%
  summarise(average_fst = sum(alpha)/sum(beta)) %>%
  ggplot(aes(x=as.factor(sample_size), y=as.factor(coverage), fill=average_fst, label=round(average_fst, 4))) +
  geom_tile() +
  geom_text() +
  scale_fill_viridis_c() +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_average_fst_raw.png", average_fst_plot, height = 4, width=6, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_average_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_average_fst_raw.png)<!-- -->

  - Although Fst is quite consistently estimated, these estimations are
    much higher than the true value (\~0.013). This may cause problems
    when estimated Fst is used to infer demography, such as migration
    rate between populations. The relative value of these genome-wide
    Fst, however, may still be trusted.

  - I also tried to get average Fst only from neutral regions on the
    genome, but that didn’t help much. At high coverage and high sample
    size, the estimated Fst is around 0.029, still much higher than the
    true value.

  - Also, at smaller sample size, Fst tends to be even more
    overestimated. This is consistent with empircal data. But
    couterintuitively, higher coverage makes the problem worse in such
    cases. We’ll need to look into this further.

### Plot the estimated per-SNP Fst (with no minimum individual filter)

``` r
fst_plot <- ggplot(fst_n_ind_final, aes(x=position, y=fst)) +
  geom_point(alpha=0.1, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_fst_raw.png", fst_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_fst_raw.png)<!-- -->

### Plot genome-wide average Fst (with minimum individual filter)

I am doing this to help Matt check whether filtering can introduce
systematic changes in genome-wide average Fst estimation. There doesn’t
appear to be any systematic bias caused by filtering in here. However,
mapping is not simulated here, and differential mapping at high Fst
regions may cause systematic biases.

``` r
fst_n_ind_final_filtered <- group_by(fst_n_ind_final, coverage, sample_size) %>%
  filter(p1_n_ind >= quantile(p1_n_ind)[4], p2_n_ind >= quantile(p2_n_ind)[4]) # filtering n_ind by the third quantile
count(fst_n_ind_final_filtered) %>%
  pivot_wider(names_from = sample_size, values_from = n)
average_fst_plot <- fst_n_ind_final_filtered %>%
  summarise(average_fst = sum(alpha)/sum(beta)) %>%
  ggplot(aes(x=as.factor(sample_size), y=as.factor(coverage), fill=average_fst, label=round(average_fst, 4))) +
  geom_tile() +
  geom_text() +
  scale_fill_viridis_c() +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_filtered_average_fst_raw.png", average_fst_plot, height = 4, width=6, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_filtered_average_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_filtered_average_fst_raw.png)<!-- -->

### Plot the estimated per-SNP Fst (with minimum individual filter)

``` r
filtered_fst_plot <- fst_n_ind_final_filtered %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.1, size=0.1) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
    facet_grid(coverage~sample_size) +
    theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_filtered_fst_raw.png", filtered_fst_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_filtered_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_filtered_fst_raw.png)<!-- -->

### Compute and plot the estimated windowed Fst (with no minimum individual filter and 1,000bp fixed windows)

``` r
windowed_fst_plot <- fixed_windowed_fst(fst_n_ind_final, 1000) %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.5, size=0.1) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
    facet_grid(coverage~sample_size) +
    theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_windowed_fst_raw.png", windowed_fst_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_windowed_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_windowed_fst_raw.png)<!-- -->

### Selection scan using PCAngsd

``` r
selection_scan <- get_selection_scan("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/")
selection_scan_summary <- group_by(selection_scan, coverage, sample_size) %>%
  summarize(n_snp = n(), log_p_cutoff = -log(0.05/n_snp))
selection_scan_plot <- ggplot(selection_scan, aes(x = pos, y = neg_log_p_value)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=27), color="red", size=0.2, shape=8) +
  geom_hline(data = selection_scan_summary, aes(yintercept = log_p_cutoff), linetype = "dashed") +
  facet_grid(coverage ~ sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_selection_scan.png", selection_scan_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_selection_scan.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_selection_scan.png)<!-- -->

## Inference with GATK’s GL model

### Fst in 1000bp windows

Note that these saf files were generated without a minMaf or a snp\_pval
filter, but a minDepth filters.

##### Read in per-SNP Fst and caculate average in 1000bp windows

``` r
fst_windowed_final_gatk <- get_estimated_windowed_fst_gatk("../two_pop_sim_fixed_m2_pos/rep_1/", 1000)
windowed_fst_plot <- fst_windowed_final_gatk %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.5, size=0.1) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
    facet_grid(coverage~sample_size) +
    theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_windowed_fst_gatk.png", windowed_fst_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_windowed_fst_gatk.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_windowed_fst_gatk.png)<!-- -->

##### Chromosome-average Fst

``` r
average_fst_gatk_plot <- fst_windowed_final_gatk %>%
  group_by(coverage, sample_size) %>%
  summarize(average_fst = sum(alpha) / sum(beta)) %>%
  ggplot(aes(x=as.factor(sample_size), y=as.factor(coverage), fill=average_fst, label=round(average_fst, 4))) +
  geom_tile() +
  geom_text() +
  scale_fill_viridis_c() +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_average_fst_gatk.png", average_fst_gatk_plot, height = 4, width=6, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_average_fst_gatk.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_average_fst_gatk.png)<!-- -->

In general, these Fst results are very similar to those from the
Samtools GL model.

### Selection scan using neutrality test stats

#### Read in the data

``` r
p1_neutrality_stats <- get_neutrality_stats("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/")
```

#### Genome wide stats

###### SFS

``` r
i=1
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40, 80, 160)){
    sfs <- scan(paste0("../two_pop_sim_fixed_m2_pos/rep_1/angsd_gatk/bam_list_p1_", sample_size, "_", coverage, "x.sfs")) %>%
      enframe(name = frequency) %>%
      mutate(frequency=(0:(sample_size*2))/(sample_size*2), coverage=coverage, sample_size=sample_size)
    if (i==1){
      sfs_final <- sfs
    } else {
      sfs_final <- bind_rows(sfs_final, sfs)
    }
    i=i+1
  }
}
sfs_final_sum <- filter(sfs_final, frequency>0, frequency<1) %>%
  group_by(coverage, sample_size) %>%
  summarise(n=sum(value))
filter(sfs_final, frequency>0, frequency<1) %>%
  group_by(coverage, sample_size) %>%
  ggplot(aes(x=frequency, y=value)) +
  geom_point(size=1) +
  geom_line() +
  geom_text(data=sfs_final_sum, x=0.8, y=40000, aes(label=paste0("n=",round(n,0)))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

###### Tajima’s estimator

``` r
p1_neutrality_stats %>%
  group_by(coverage, sample_size) %>%
  summarise(t_p = round(sum(t_p) / sum(n_sites), 5)) %>%
  pivot_wider(names_from = sample_size, values_from = t_p) %>% 
  kable()
```

| coverage |       5 |      10 |      20 |      40 |      80 |     160 |
| -------: | ------: | ------: | ------: | ------: | ------: | ------: |
|     0.25 | 0.00393 | 0.00382 | 0.00385 | 0.00383 | 0.00385 | 0.00387 |
|     0.50 | 0.00388 | 0.00381 | 0.00385 | 0.00384 | 0.00386 | 0.00389 |
|     1.00 | 0.00387 | 0.00381 | 0.00386 | 0.00386 | 0.00387 | 0.00388 |
|     2.00 | 0.00385 | 0.00379 | 0.00386 | 0.00384 | 0.00386 | 0.00387 |
|     4.00 | 0.00385 | 0.00379 | 0.00385 | 0.00384 | 0.00384 | 0.00384 |
|     8.00 | 0.00383 | 0.00377 | 0.00381 | 0.00379 | 0.00378 | 0.00380 |

###### Watterson’s estimator

``` r
p1_neutrality_stats %>%
  group_by(coverage, sample_size) %>%
  summarise(t_w = round(sum(t_w) / sum(n_sites), 5)) %>%
  pivot_wider(names_from = sample_size, values_from = t_w) %>% 
  kable()
```

| coverage |       5 |      10 |      20 |      40 |      80 |     160 |
| -------: | ------: | ------: | ------: | ------: | ------: | ------: |
|     0.25 | 0.00414 | 0.00339 | 0.00415 | 0.00483 | 0.00596 | 0.00839 |
|     0.50 | 0.00386 | 0.00384 | 0.00426 | 0.00486 | 0.00618 | 0.00888 |
|     1.00 | 0.00391 | 0.00398 | 0.00443 | 0.00506 | 0.00641 | 0.00894 |
|     2.00 | 0.00399 | 0.00406 | 0.00456 | 0.00520 | 0.00647 | 0.00891 |
|     4.00 | 0.00403 | 0.00411 | 0.00461 | 0.00520 | 0.00641 | 0.00876 |
|     8.00 | 0.00401 | 0.00407 | 0.00454 | 0.00506 | 0.00612 | 0.00823 |

It is strange that Watterson’s estimator is now overestimated.

#### Window-based stats

#### Tajima’s D

``` r
p1_tajima_plot <- ggplot(p1_neutrality_stats, aes(x = pos, y = tajima)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=1.5), color="red", size=0.2, shape=8) +
  facet_grid(coverage ~ sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_p1_tajima_d.png", p1_tajima_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_p1_tajima_d.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_p1_tajima_d.png)<!-- -->

#### Fay and Wu’s H

``` r
max_fayh <- group_by(p1_neutrality_stats, coverage, sample_size) %>%
  summarise(max_fayh = max(fayh)) %>%
  crossing(mutations_final_m2[,2])

p1_fayh_plot <- ggplot(p1_neutrality_stats, aes(x = pos, y = fayh)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_point(data=max_fayh, aes(x=position, y=max_fayh*1.1), color="red", size=0.2, shape=8) +
  facet_wrap(coverage ~ sample_size, scales = "free_y") +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_p1_fay_h.png", p1_fayh_plot, height = 10, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_p1_fay_h.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_p1_fay_h.png)<!-- -->

## RAD seq simulaiton and inference

To simulate RAD-seq, I assumed that the genotype calling is perfectly
accurate (i.e. sequence depth is high). I then took random samples along
the chromosome representing RAD tags.

### Get true sample allele count

``` r
get_sample_allele_count_per_pop("../two_pop_sim_fixed_m2_pos/rep_1/")
fst_all_snps <- allele_count_to_fst("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/")
write_tsv(fst_all_snps, "/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/fst_all_snps.tsv")
```

### Get the true SFS and theta estimators with population samples in p1

``` r
fst_all_snps <- read_tsv("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos/rep_1/slim/fst_all_snps.tsv")
## SFS
group_by(fst_all_snps, sample_size) %>%
  filter(p1>0, p1<1) %>%
  count(p1) %>%
  ggplot(aes(x=p1, y=n)) +
  geom_line() +
  geom_point() +
  theme_cowplot() +
  facet_wrap(~sample_size)
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
## theta_p
group_by(fst_all_snps, sample_size) %>%
  filter(p1>0, p1<1) %>%
  summarise(theta_p = 2 * sum(p1*(1-p1)) / 30000000,
            theta_w = n() / 30000000 / sum(1/1:(2*sample_size-1)))
```

    ## # A tibble: 6 x 3
    ##   sample_size theta_p theta_w
    ##         <dbl>   <dbl>   <dbl>
    ## 1           5 0.00325 0.00361
    ## 2          10 0.00342 0.00359
    ## 3          20 0.00352 0.00355
    ## 4          40 0.00356 0.00351
    ## 5          80 0.00358 0.00344
    ## 6         160 0.00359 0.00336

### Get sample true MAF and Fst from allele counts

Each RAD tag is 150 bp. Note that `n_rad_tag` are numbers of RAD tags
per Mbp. According to the “Breaking RAD” paper, the median RAD tag
density that they’ve found in studies published by then is 4.08. A few
studies had up to 20 tags per Mbp, three had up to 110 per Mbp, and one
had 362 tags per Mbp.

``` r
i <- 1
for (n in c(1,2,4,8,16,32)*120){
  maf <- subset_fst(n_rad_tag = n, combined_fst = fst_all_snps) %>%
    mutate(n_rad_tag = n/30)
  if (i == 1){
    maf_final <- maf
  } else {
    maf_final <- bind_rows(maf_final, maf)
  }
  i <- i + 1
}
```

### Chromosome-wide stats

#### Plot SFS in p1 from RAD data

``` r
count(maf_final, p1, sample_size, n_rad_tag) %>%
  filter(p1>0, p1<1) %>%
  ggplot(aes(x=p1, y=n)) +
  geom_point(size=0.5) +
  geom_line(size=0.2) +
  facet_grid(n_rad_tag~sample_size, scales = "free_y") +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

#### Tajima’s estimator

``` r
group_by(maf_final, sample_size, n_rad_tag) %>%
  summarise(theta_t = 2 * sum(p1 * (1 - p1)) / unique(n_sites)) %>%
  pivot_wider(names_from = sample_size, values_from = theta_t) %>%
  kable()
```

    ## `summarise()` regrouping output by 'sample_size' (override with `.groups` argument)

| n\_rad\_tag |         5 |        10 |        20 |        40 |        80 |       160 |
| ----------: | --------: | --------: | --------: | --------: | --------: | --------: |
|           4 | 0.0033753 | 0.0035577 | 0.0036674 | 0.0037537 | 0.0038450 | 0.0038126 |
|           8 | 0.0034398 | 0.0035614 | 0.0036205 | 0.0036860 | 0.0037190 | 0.0037028 |
|          16 | 0.0033141 | 0.0034242 | 0.0035067 | 0.0035776 | 0.0036231 | 0.0036235 |
|          32 | 0.0032981 | 0.0034204 | 0.0035125 | 0.0035752 | 0.0036206 | 0.0036185 |
|          64 | 0.0032879 | 0.0034148 | 0.0035083 | 0.0035587 | 0.0035994 | 0.0036082 |
|         128 | 0.0032864 | 0.0034177 | 0.0035249 | 0.0035683 | 0.0036006 | 0.0036114 |

#### Watterson’s estimator

``` r
filter(maf_final, p1 > 0, p1 < 1) %>%
  group_by(sample_size, n_rad_tag) %>%
  summarize(theta_w = n() / (unique(n_sites)*sum(1/(1:(unique(sample_size)*2-1))))) %>%
  pivot_wider(names_from = sample_size, values_from = theta_w) %>%
  kable()
```

    ## `summarise()` regrouping output by 'sample_size' (override with `.groups` argument)

| n\_rad\_tag |         5 |        10 |        20 |        40 |        80 |       160 |
| ----------: | --------: | --------: | --------: | --------: | --------: | --------: |
|           4 | 0.0038236 | 0.0039045 | 0.0038275 | 0.0037550 | 0.0036634 | 0.0035928 |
|           8 | 0.0038138 | 0.0037723 | 0.0036588 | 0.0035934 | 0.0034973 | 0.0034449 |
|          16 | 0.0036377 | 0.0036093 | 0.0035365 | 0.0034945 | 0.0034428 | 0.0033750 |
|          32 | 0.0036604 | 0.0036004 | 0.0035632 | 0.0034952 | 0.0034374 | 0.0033648 |
|          64 | 0.0036750 | 0.0036155 | 0.0035690 | 0.0035026 | 0.0034251 | 0.0033638 |
|         128 | 0.0036608 | 0.0035999 | 0.0035741 | 0.0035192 | 0.0034381 | 0.0033565 |

### Plot per SNP Fst

``` r
mutate(maf_final, coverage="RAD") %>%
  filter(maf_mean>0.05, maf_mean < 0.95) %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.5, size=0.2) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
    facet_grid(n_rad_tag~sample_size) +
    theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

# Two populations with divergent selection, with smaller population size ( Ne\~10,000 in each population)

The same population size is simulated, but I’ve scaled down mutation
rate, recombination rate, migration rate. The selection coefficient is
unchanged. (Ignore the directory name. I’m just using it for
convenience.)

## The model

### Read in the ancestral states

``` r
ancestral <- get_ancestral("../two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
```

### Read mutation and substitution file

``` r
mutations_final <- get_mutations("../two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
real_theta_t_p1 <- sum(2*mutations_final$p1*(1-mutations_final$p1))/30000000
real_theta_t_p2 <- sum(2*mutations_final$p2*(1-mutations_final$p2))/30000000
real_theta_w_p1 <- filter(mutations_final, p1 > 0, p1 < 1) %>%
  nrow() %>% `/`(30000000*sum(1/(1:(5000*2-1))))
real_theta_w_p2 <- filter(mutations_final, p2 > 0, p2 < 1) %>%
  nrow() %>% `/`(30000000*sum(1/(1:(5000*2-1))))
tibble(theta = c("tajima", "watterson"), p1=c(real_theta_t_p1, real_theta_w_p1), p2=c(real_theta_t_p2,real_theta_w_p2))
```

    ## # A tibble: 2 x 3
    ##   theta           p1       p2
    ##   <chr>        <dbl>    <dbl>
    ## 1 tajima    0.000704 0.000728
    ## 2 watterson 0.000593 0.000605

### Plot Fst

``` r
mutations_final_m1 <- filter(mutations_final, type=="m1")
mutations_final_m2 <- filter(mutations_final, type=="m2")
mean_neutral_fst_weighted <- sum(mutations_final_m1$h_t-mutations_final_m1$h_s)/sum(mutations_final_m1$h_t)
mean_neutral_fst_weighted
```

    ## [1] 0.03674846

``` r
ggplot(mutations_final_m1, aes(x=position/10^6, y=fst)) +
  geom_point(size=0.02, alpha=0.5) +
  geom_point(data=mutations_final_m2, color="red") +
  labs(x = "position (in Mbp)", y = expression(F[ST])) + 
  xlim(c(0, 30)) +
  ylim(c(0, 1)) +
  theme_cowplot() +
  theme(text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

``` r
arrange(mutations_final, desc(fst)) %>%
  head()
```

    ## # A tibble: 6 x 9
    ##   type  position ancestral    p1     p2 frequency_mean   h_t    h_s   fst
    ##   <chr>    <dbl> <chr>     <dbl>  <dbl>          <dbl> <dbl>  <dbl> <dbl>
    ## 1 m2     2500001 G         0.969 0.0019          0.485 0.500 0.0320 0.936
    ## 2 m1     2498645 C         0.970 0.0218          0.496 0.500 0.0506 0.899
    ## 3 m1     2491161 A         0.933 0.0326          0.483 0.499 0.0940 0.812
    ## 4 m1     2477186 A         0.887 0.0084          0.448 0.495 0.108  0.781
    ## 5 m1     2469162 A         0.879 0.0084          0.444 0.494 0.115  0.767
    ## 6 m1     2470141 A         0.879 0.0084          0.444 0.494 0.115  0.767

``` r
arrange(mutations_final_m2, desc(fst))
```

    ## # A tibble: 8 x 9
    ##   type  position ancestral    p1     p2 frequency_mean   h_t    h_s   fst
    ##   <chr>    <dbl> <chr>     <dbl>  <dbl>          <dbl> <dbl>  <dbl> <dbl>
    ## 1 m2     2500001 G         0.969 0.0019          0.485 0.500 0.0320 0.936
    ## 2 m2    22500001 G         0.869 0.0082          0.439 0.492 0.122  0.752
    ## 3 m2    17500001 A         0.833 0.002           0.417 0.486 0.141  0.709
    ## 4 m2     7500001 A         0.766 0.0018          0.384 0.473 0.181  0.618
    ## 5 m2     5000001 G         0.713 0.0009          0.357 0.459 0.205  0.553
    ## 6 m2    15000001 C         0.645 0.002           0.324 0.438 0.231  0.472
    ## 7 m2    25000001 G         0.582 0.003           0.292 0.414 0.246  0.404
    ## 8 m2    12500001 A         0.299 0.0008          0.150 0.255 0.211  0.175

## Inference based on Samtool’s GL model

### Read in read depth and estimated Fst

``` r
fst_n_ind_final <- get_estimated_fst("../two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
```

### Plot the estimated per-SNP Fst (with no minimum individual filter)

``` r
fst_plot <- ggplot(fst_n_ind_final, aes(x=position, y=fst)) +
  geom_point(alpha=0.1, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=0.2, shape=8) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_fst_raw.png", fst_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_fst_raw.png)<!-- -->

### Compute and plot the estimated windowed Fst (with no minimum individual filter and 5,000bp fixed windows)

``` r
windowed_fst_plot <- fixed_windowed_fst(fst_n_ind_final, 5000) %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.8, size=0.15) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=2, shape=8) +
    ylim(NA, 1.05) +
    facet_grid(coverage~sample_size) +
    theme_cowplot() +
    theme(text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_windowed_fst_raw.png", windowed_fst_plot, width = 42, height = 21, units = "cm", pointsize = 20)
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_windowed_fst_raw.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_windowed_fst_raw.png)<!-- -->

### Selection scan using PCAngsd

``` r
selection_scan <- get_selection_scan("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
selection_scan_summary <- group_by(selection_scan, coverage, sample_size) %>%
  summarize(n_snp = n(), log_p_cutoff = -log(0.05/n_snp))
selection_scan_plot <- ggplot(selection_scan, aes(x = pos, y = neg_log_p_value)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=16), color="red", size=0.2, shape=8) +
  geom_hline(data = selection_scan_summary, aes(yintercept = log_p_cutoff), linetype = "dashed") +
  facet_grid(coverage ~ sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection_scan.png", selection_scan_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection_scan.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_selection_scan.png)<!-- -->

## Inference based on GATK’s GL model

### Selection scan using neutrality test stats

#### Read in the data

``` r
p1_neutrality_stats <- get_neutrality_stats("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
```

#### Genome wide stats

###### SFS

``` r
i=1
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40, 80, 160)){
    sfs <- scan(paste0("../two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/angsd_gatk/bam_list_p1_", sample_size, "_", coverage, "x.sfs")) %>%
      enframe(name = frequency) %>%
      mutate(frequency=(0:(sample_size*2))/(sample_size*2), coverage=coverage, sample_size=sample_size)
    if (i==1){
      sfs_final <- sfs
    } else {
      sfs_final <- bind_rows(sfs_final, sfs)
    }
    i=i+1
  }
}
sfs_final_sum <- filter(sfs_final, frequency>0, frequency<1) %>%
  group_by(coverage, sample_size) %>%
  summarise(n=sum(value))
filter(sfs_final, frequency>0, frequency<1) %>%
  group_by(coverage, sample_size) %>%
  ggplot(aes(x=frequency, y=value)) +
  geom_point(size=0.5) +
  geom_line() +
  geom_text(data=sfs_final_sum, x=0.8, y=40000, aes(label=paste0("n=",round(n,0)))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

###### Tajima’s estimator

``` r
p1_neutrality_stats %>%
  group_by(coverage, sample_size) %>%
  summarise(t_p = round(sum(t_p) / sum(n_sites), 5)) %>%
  pivot_wider(names_from = sample_size, values_from = t_p)
```

    ## # A tibble: 6 x 7
    ## # Groups:   coverage [6]
    ##   coverage     `5`    `10`    `20`    `40`    `80`   `160`
    ##      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1     0.25 0.00106 0.00118 0.00124 0.00116 0.00113 0.00112
    ## 2     0.5  0.00119 0.00134 0.00125 0.00117 0.00115 0.00114
    ## 3     1    0.00135 0.00138 0.00126 0.00120 0.00117 0.00117
    ## 4     2    0.00152 0.00138 0.00134 0.00126 0.00123 0.00122
    ## 5     4    0.00148 0.00148 0.00143 0.00134 0.0013  0.00129
    ## 6     8    0.0015  0.00153 0.00147 0.00136 0.00132 0.00131

###### Watterson’s estimator

``` r
p1_neutrality_stats %>%
  group_by(coverage, sample_size) %>%
  summarise(t_w = round(sum(t_w) / sum(n_sites), 5)) %>%
  pivot_wider(names_from = sample_size, values_from = t_w)
```

    ## # A tibble: 6 x 7
    ## # Groups:   coverage [6]
    ##   coverage     `5`    `10`    `20`    `40`    `80`  `160`
    ##      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ## 1     0.25 0.00087 0.00104 0.00181 0.00329 0.00568 0.0100
    ## 2     0.5  0.00112 0.00108 0.00191 0.00389 0.00627 0.0107
    ## 3     1    0.00124 0.00161 0.00304 0.00438 0.0069  0.0116
    ## 4     2    0.00142 0.00248 0.00352 0.00497 0.00778 0.0130
    ## 5     4    0.00205 0.00285 0.00404 0.00568 0.00888 0.0149
    ## 6     8    0.00211 0.003   0.00427 0.00596 0.00930 0.0157

Again, Watterson’s estimator is now severely overestimated at higher
sample size and higher coverage.

#### Window-based stats

###### Plot Tajima’s D

``` r
p1_tajima_plot <- ggplot(p1_neutrality_stats, aes(x = pos, y = tajima)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_point(data=mutations_final_m2, aes(x=position, y=1.5), color="red", size=0.2, shape=8) +
  facet_grid(coverage ~ sample_size) +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_tajima_d.png", p1_tajima_plot, height = 8, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_tajima_d.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_tajima_d.png)<!-- -->

###### Plot Fay and Wu’s H

``` r
max_fayh <- group_by(p1_neutrality_stats, coverage, sample_size) %>%
  summarise(max_fayh = max(fayh)) %>%
  crossing(mutations_final_m2[,2])

p1_fayh_plot <- ggplot(p1_neutrality_stats, aes(x = pos, y = fayh)) +
  geom_point(alpha=0.5, size=0.1) +
  geom_smooth(span = 0.01) +
  geom_point(data=max_fayh, aes(x=position, y=max_fayh*1.1), color="red", size=0.2, shape=8) +
  facet_wrap(coverage ~ sample_size, scales = "free_y") +
  theme_cowplot()
ggsave("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_fay_h.png", p1_fayh_plot, height = 10, width=15, units = "in")
```

``` r
include_graphics("../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_fay_h.png")
```

![](../figures/two_pop_sim_fixed_m2_pos_lower_s_lower_r_p1_fay_h.png)<!-- -->

## RAD seq simulation and inference

To simulate RAD-seq, I assumed that the genotype calling is perfectly
accurate (i.e. sequence depth is high). I then took random samples along
the chromosome representing RAD tags.

### Get true sample allele count

``` r
get_sample_allele_count_per_pop("../two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
fst_all_snps <- allele_count_to_fst("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/")
write_tsv(fst_all_snps, "/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/slim/fst_all_snps.tsv")
```

### Get the true SFS and theta estimators with population samples in p1

``` r
fst_all_snps <- read_tsv("/workdir/lcwgs-simulation/two_pop_sim_fixed_m2_pos_lower_s_lower_r/rep_1/slim/fst_all_snps.tsv")
## SFS
group_by(fst_all_snps, sample_size) %>%
  filter(p1>0, p1<1) %>%
  count(p1) %>%
  ggplot(aes(x=p1, y=n)) +
  geom_line() +
  geom_point() +
  theme_cowplot() +
  facet_wrap(~sample_size)
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-58-1.png)<!-- -->

``` r
## theta_p
group_by(fst_all_snps, sample_size) %>%
  filter(p1>0, p1<1) %>%
  summarise(theta_p = 2 * sum(p1*(1-p1)) / 30000000,
            theta_w = n() / 30000000 / sum(1/1:(2*sample_size-1)))
```

    ## # A tibble: 6 x 3
    ##   sample_size  theta_p  theta_w
    ##         <dbl>    <dbl>    <dbl>
    ## 1           5 0.000630 0.000693
    ## 2          10 0.000666 0.000688
    ## 3          20 0.000685 0.000683
    ## 4          40 0.000694 0.000669
    ## 5          80 0.000698 0.000652
    ## 6         160 0.000701 0.000633

### Get sample true MAF and Fst from allele counts

Each RAD tag is 150 bp. Note that `n_rad_tag` are numbers of RAD tags
per Mbp. According to the “Breaking RAD” paper, the median RAD tag
density that they’ve found in studies published by then is 4.08. A few
studies had up to 20 tags per Mbp, three had up to 110 per Mbp, and one
had 362 tags per Mbp.

``` r
i <- 1
for (n in c(1,2,4,8,16,32)*120){
  maf <- subset_fst(n_rad_tag = n, combined_fst = fst_all_snps) %>%
    mutate(n_rad_tag = n/30)
  if (i == 1){
    maf_final <- maf
  } else {
    maf_final <- bind_rows(maf_final, maf)
  }
  i <- i + 1
}
```

### Plot SFS in p1 from RAD data

``` r
count(maf_final, p1, sample_size, n_rad_tag) %>%
  filter(p1>0, p1<1) %>%
  ggplot(aes(x=p1, y=n)) +
  geom_point(size=0.5) +
  geom_line(size=0.2) +
  facet_grid(n_rad_tag~sample_size, scales = "free_y") +
  theme_cowplot()
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

### Get theta estimation in p1 from RAD data

#### Tajima’s estimator

``` r
group_by(maf_final, sample_size, n_rad_tag) %>%
  summarise(theta_t = 2 * sum(p1 * (1 - p1)) / unique(n_sites)) %>%
  pivot_wider(names_from = sample_size, values_from = theta_t) %>%
  kable()
```

    ## `summarise()` regrouping output by 'sample_size' (override with `.groups` argument)

| n\_rad\_tag |         5 |        10 |        20 |        40 |        80 |       160 |
| ----------: | --------: | --------: | --------: | --------: | --------: | --------: |
|           4 | 0.0006854 | 0.0007307 | 0.0007570 | 0.0007569 | 0.0007707 | 0.0007692 |
|           8 | 0.0006932 | 0.0007463 | 0.0007755 | 0.0007788 | 0.0007943 | 0.0007970 |
|          16 | 0.0006821 | 0.0007351 | 0.0007532 | 0.0007536 | 0.0007611 | 0.0007696 |
|          32 | 0.0006350 | 0.0006794 | 0.0006992 | 0.0007023 | 0.0007053 | 0.0007118 |
|          64 | 0.0006380 | 0.0006707 | 0.0006789 | 0.0006886 | 0.0006918 | 0.0006972 |
|         128 | 0.0006342 | 0.0006650 | 0.0006746 | 0.0006857 | 0.0006865 | 0.0006893 |

#### Watterson’s estimator

``` r
filter(maf_final, p1 > 0, p1 < 1) %>%
  group_by(sample_size, n_rad_tag) %>%
  summarize(theta_w = n() / (unique(n_sites)*sum(1/(1:(unique(sample_size)*2-1))))) %>%
  pivot_wider(names_from = sample_size, values_from = theta_w) %>%
  kable()
```

    ## `summarise()` regrouping output by 'sample_size' (override with `.groups` argument)

| n\_rad\_tag |         5 |        10 |        20 |        40 |        80 |       160 |
| ----------: | --------: | --------: | --------: | --------: | --------: | --------: |
|           4 | 0.0007218 | 0.0007311 | 0.0008044 | 0.0007577 | 0.0007815 | 0.0007394 |
|           8 | 0.0007608 | 0.0007389 | 0.0007785 | 0.0007354 | 0.0007180 | 0.0006872 |
|          16 | 0.0007520 | 0.0007281 | 0.0007534 | 0.0007279 | 0.0006969 | 0.0006924 |
|          32 | 0.0006911 | 0.0006913 | 0.0007146 | 0.0007099 | 0.0006848 | 0.0006719 |
|          64 | 0.0006988 | 0.0006967 | 0.0006861 | 0.0006717 | 0.0006514 | 0.0006325 |
|         128 | 0.0006884 | 0.0006833 | 0.0006824 | 0.0006655 | 0.0006438 | 0.0006205 |

### Plot per SNP Fst

``` r
mutate(maf_final, coverage="RAD") %>%
  filter(maf_mean>0.05, maf_mean < 0.95) %>%
  ggplot(aes(x=position, y=fst)) +
    geom_point(alpha=0.8, size=0.15) +
    geom_point(data=mutations_final_m2, aes(x=position, y=1.01), color="red", size=2, shape=8) +
    ylim(NA, 1.05) +
    facet_grid(n_rad_tag~sample_size) +
    theme_cowplot() +
    theme(text = element_text(size=20),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
```

![](data_analysis_two_pop_fixed_m2_pos_files/figure-gfm/unnamed-chunk-63-1.png)<!-- -->
