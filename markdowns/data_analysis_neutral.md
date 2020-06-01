Data analysis with neutral simulation
================

  - [Define some functions](#define-some-functions)
  - [Data wrangling with SLiM output](#data-wrangling-with-slim-output)
      - [Read in the ancestral states](#read-in-the-ancestral-states)
      - [Read mutation and substitution
        file](#read-mutation-and-substitution-file)
      - [Data wrangling with the mutation
        file](#data-wrangling-with-the-mutation-file)
      - [Plot the true site frequency
        spectrum](#plot-the-true-site-frequency-spectrum)
  - [ANGSD results](#angsd-results)
      - [Read maf estimation and join mutation and maf
        files](#read-maf-estimation-and-join-mutation-and-maf-files)
      - [Plot the estimated site frequency
        spectrum](#plot-the-estimated-site-frequency-spectrum)
      - [Plot the estimated allele frequency
        distribution](#plot-the-estimated-allele-frequency-distribution)
      - [Plot estimated allele frequency vs. true allele frequency (this
        includes the false positives but not the false
        negatives)](#plot-estimated-allele-frequency-vs.-true-allele-frequency-this-includes-the-false-positives-but-not-the-false-negatives)
      - [Plot estimated allele frequency vs. true allele frequency in
        bins (this includes the false positives but not the false
        negatives)](#plot-estimated-allele-frequency-vs.-true-allele-frequency-in-bins-this-includes-the-false-positives-but-not-the-false-negatives)
      - [Plot error vs. true allele frequency in
        bins](#plot-error-vs.-true-allele-frequency-in-bins)
      - [Check the SNPs with highest
        error](#check-the-snps-with-highest-error)
      - [True frequency distribution of false
        negatives](#true-frequency-distribution-of-false-negatives)
      - [Esimated frequency distribution of false
        positives](#esimated-frequency-distribution-of-false-positives)
      - [Read windowed thetas estimated by
        `realSFS`](#read-windowed-thetas-estimated-by-realsfs)
      - [real theta values with the entire
        population](#real-theta-values-with-the-entire-population)
      - [thetas estimated from the estimated allele frequencies (not
        from
        `realSFS`)](#thetas-estimated-from-the-estimated-allele-frequencies-not-from-realsfs)
      - [thetas estimated from
        `realSFS`](#thetas-estimated-from-realsfs)
      - [Plot Watterson’s estimator and Tajima’s estimator of theta and
        Tajima’s D in 10,000bp fixed
        windows](#plot-wattersons-estimator-and-tajimas-estimator-of-theta-and-tajimas-d-in-10000bp-fixed-windows)
  - [Compare individual barcoding with
    Pool-seq](#compare-individual-barcoding-with-pool-seq)
      - [Get minor allele frequencies estimated from
        Pool-seq](#get-minor-allele-frequencies-estimated-from-pool-seq)
      - [Plot estimated allele frequency vs. true allele frequency (this
        includes the false positives but not the false
        negatives)](#plot-estimated-allele-frequency-vs.-true-allele-frequency-this-includes-the-false-positives-but-not-the-false-negatives-1)
      - [Plot estimated allele frequency vs. true allele frequency in
        bins (this includes the false positives but not the false
        negatives)](#plot-estimated-allele-frequency-vs.-true-allele-frequency-in-bins-this-includes-the-false-positives-but-not-the-false-negatives-1)
      - [Plot absolute values of error vs. true allele frequency in bins
        (this includes the false positives but not the false
        negatives)](#plot-absolute-values-of-error-vs.-true-allele-frequency-in-bins-this-includes-the-false-positives-but-not-the-false-negatives)

``` r
library(tidyverse)
library(cowplot)
library(knitr)
```

# Define some functions

``` r
summarise_by_design <- function(joined_frequency_table){
  joined_frequency_table %>%
  mutate(error_squared=error^2) %>%
  group_by(coverage, sample_size) %>%
  summarise(r_squared=paste0("R^2==", round(summary(lm(estimated_frequency~frequency))$r.squared,3)), 
            n=paste0("n==",n()), 
            root_mean_error_squared=paste0("RMSE==", round(sqrt(mean(error_squared)), 3)))
}
plot_frequency <- function(joined_frequency_table, joined_summary_table){
  joined_frequency_table %>%
    ggplot(aes(x=frequency, y=estimated_frequency)) +
    geom_point(alpha=0.1, size=0.1) +
    geom_smooth(method="lm", color="red", size=1, se = F) +
    geom_text(data = joined_summary_table, x = 0.86, y = 0.25, aes(label=r_squared), color = 'red',  parse = TRUE) +
    geom_text(data = joined_summary_table, x = 0.86, y = 0.12, aes(label=n), color = 'red',  parse = TRUE) +
    facet_grid(coverage~sample_size) +
    theme_cowplot()
}
plot_frequency_in_bins <- function(joined_frequency_table, joined_summary_table){
  joined_frequency_table %>%
    ggplot(aes(x=frequency_bin, y=estimated_frequency)) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(data = joined_summary_table, x = 8.6, y = 0.2, aes(label=r_squared), color = 'black',  parse = TRUE) +
    geom_text(data = joined_summary_table, x = 8.6, y = 0.07, aes(label=n), color = 'black',  parse = TRUE) +
    scale_x_discrete(labels=seq(0.05, 0.95, 0.1))  +
    facet_grid(coverage~sample_size) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45))
}
plot_error_in_bins <- function(joined_frequency_table, joined_summary_table){
  joined_frequency_table %>%
    ggplot(aes(x=frequency_bin, y=abs(error))) +
    geom_boxplot(outlier.shape = NA) +
    geom_text(data = joined_summary_table, x = 8, y = 0.77, aes(label=root_mean_error_squared), color = 'black',  parse = TRUE) +
    geom_text(data = joined_summary_table, x = 8, y = 0.64, aes(label=n), color = 'black',  parse = TRUE) +
    facet_grid(coverage~sample_size) +
    scale_x_discrete(labels=seq(0.05, 0.95, 0.1))  +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45))
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
```

# Data wrangling with SLiM output

## Read in the ancestral states

``` r
ancestral <- read_csv("../neutral_sim/rep_1/slim/ancestral.fasta")[[1]] %>%
  str_split(pattern="") %>%
  .[[1]] %>%
  bind_cols(ancestral=., position=1:30000000)
```

## Read mutation and substitution file

``` r
## Read in the mutation file outputted by SLiM
mutations <- read_delim("../neutral_sim/rep_1/slim/mutations.txt", delim = " ", col_names = F) %>%
  transmute(type=X6, position=X7+1, base=X13, frequency=X12/2000) %>%
  left_join(ancestral, by="position") %>%
  group_by(type, position, ancestral, base) %>%
  summarise(frequency=sum(frequency)) %>%
  ungroup()
## Read in the substitutions file outputted by SLiM
## This is necessary because mutations can happen again after one fixation, so frequencies from the mutation file do not always reflect the true derived allele frequency
substitutions <- read_delim("../neutral_sim/rep_1/slim/substitutions.txt", delim = " ", skip=2, col_names = F) %>%
  transmute(type=X3, position=X4+1, base=X10, generation=X9) %>%
  group_by(type, position) %>%
  filter(generation==max(generation)) %>%
  ungroup() %>%
  left_join(ancestral, by="position") %>%
  select(-generation) %>%
  filter(base!=ancestral) %>%
  arrange(position)
```

## Data wrangling with the mutation file

The following steps are necessary because there are complications such
as back mutations and triallelic loci in the mutation file

``` r
## Join mutations and substitutions in a temp table
mutations_final_temp_1 <-  mutations %>%
  spread(key = base, value=frequency) %>%
  full_join(substitutions, by=c("position", "type", "ancestral")) %>%
  arrange(position) %>%
  mutate(base=ifelse(is.na(base), ancestral, base)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(frequency=1-`A` -`C` -`G` -`T`)
## More wrangling
mutations_final_temp_2 <- mutations_final_temp_1[1:7] %>%
  gather(key=base, value=frequency, 4:7) %>%
  bind_rows(mutations_final_temp_1[c(1:3, 8:9)]) %>%
  mutate(frequency=ifelse(base==ancestral, 0, frequency)) %>%
  group_by(type, position, ancestral) %>%
  filter(frequency!=0) %>%
  summarise(frequency=sum(frequency), base=paste0(base, collapse = "")) %>%
  ungroup() 
mutations_final <- mutations_final_temp_2 %>%
  filter(frequency!=1)
substitutions_pos <- mutations_final_temp_2 %>%
  filter(frequency==1) %>%
  .$position %>%
  c(., substitutions$position) %>%
  unique()
```

## Plot the true site frequency spectrum

``` r
n_mutations <- dim(mutations_final)[1]
ggplot(mutations_final, aes(x=frequency)) +
  geom_histogram(bins=201) +
  annotate("text", x=0.8, y=150000, label=paste0("n=", n_mutations)) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# ANGSD results

## Read maf estimation and join mutation and maf files

``` r
i=1
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40, 80, 160)){
    ## read in estimated maf
    maf <- read_tsv(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.mafs.gz")) %>%
      mutate(estimated_frequency=knownEM) %>%
      select(position, major, minor, anc, estimated_frequency, nInd) %>%
      arrange(position)
    ## join estimated maf with true snps and only keep the snps that exist in estimated maf
    joined_frequency <- right_join(mutations_final, maf, by="position") %>%
      select(-ancestral) %>%
      mutate(coverage=coverage, sample_size=sample_size, 
             frequency=ifelse(is.na(frequency) & !(position %in% substitutions_pos), 0, frequency),
             frequency=ifelse(is.na(frequency) & position %in% substitutions_pos, 1, frequency))
    ## find false negatives
    false_negatives <- anti_join(mutations_final, maf, by="position") %>%
      mutate(coverage=coverage, sample_size=sample_size)
    ## find false positives
    false_positives <- anti_join(maf, mutations_final, by="position") %>%
      mutate(coverage=coverage, sample_size=sample_size)
    ## read in estimated sfs
    sfs <- scan(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.sfs")) %>%
      enframe(name = frequency) %>%
      mutate(frequency=(0:(sample_size*2))/(sample_size*2), coverage=coverage, sample_size=sample_size)
    ## compile the final files for plotting
    if (i==1){
      joined_frequency_final <- joined_frequency
      false_negatives_final <- false_negatives
      false_positives_final <- false_positives
      sfs_final <- sfs
    } else {
      joined_frequency_final <- bind_rows(joined_frequency_final, joined_frequency)
      false_negatives_final <- bind_rows(false_negatives_final, false_negatives)
      false_positives_final <- bind_rows(false_positives_final, false_positives)
      sfs_final <- bind_rows(sfs_final, sfs)
    }
    i=i+1
  }
}
write_tsv(joined_frequency_final, "../neutral_sim/rep_1/angsd/joined_frequency_final.tsv")
write_tsv(false_negatives_final, "../neutral_sim/rep_1/angsd/false_negatives_final.tsv")
write_tsv(false_positives_final, "../neutral_sim/rep_1/angsd/false_positives_final.tsv")
write_tsv(sfs_final, "../neutral_sim/rep_1/angsd/sfs_final.tsv")
```

``` r
joined_frequency_final <- read_tsv("../neutral_sim/rep_1/angsd/joined_frequency_final.tsv") %>%
  mutate(frequency_bin = cut(frequency, breaks = 0:10/10, include.lowest = T), error=estimated_frequency-frequency)
false_negatives_final <- read_tsv("../neutral_sim/rep_1/angsd/false_negatives_final.tsv")
false_positives_final <- read_tsv("../neutral_sim/rep_1/angsd/false_positives_final.tsv")
sfs_final <- read_tsv("../neutral_sim/rep_1/angsd/sfs_final.tsv")
```

## Plot the estimated site frequency spectrum

These are obtained from `realSFS`.

``` r
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

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
filter(sfs_final, frequency>0, frequency<1) %>%
  ggplot(aes(x=frequency, y=value)) +
  geom_col(aes(width = 0.5/sample_size)) +
  geom_text(data=sfs_final_sum, x=0.8, y=40000, aes(label=paste0("n=",round(n,0)))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

## Plot the estimated allele frequency distribution

These are the histogram of estimated allele frequencies

``` r
joined_summary <- summarise_by_design(joined_frequency_final)

joined_frequency_final %>%
  ggplot(aes(x=estimated_frequency)) +
  geom_histogram() +
  geom_text(data=joined_summary, x=0.8, y=20000, aes(label=n), parse=T) +
  facet_grid(coverage~sample_size, scales ="free_y") +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Plot estimated allele frequency vs. true allele frequency (this includes the false positives but not the false negatives)

``` r
estimated_vs_true_frequency <- plot_frequency(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_estimated_vs_true_frequency.png", estimated_vs_true_frequency, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_estimated_vs_true_frequency.png")
```

![](../figures/neutral_estimated_vs_true_frequency.png)<!-- -->

## Plot estimated allele frequency vs. true allele frequency in bins (this includes the false positives but not the false negatives)

``` r
estimated_vs_true_frequency_bin <- plot_frequency_in_bins(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_estimated_vs_true_frequency_bin.png", estimated_vs_true_frequency_bin, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_estimated_vs_true_frequency_bin.png")
```

![](../figures/neutral_estimated_vs_true_frequency_bin.png)<!-- -->

## Plot error vs. true allele frequency in bins

``` r
error_vs_true_frequency_bin <- plot_error_in_bins(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_error_vs_true_frequency_bin.png", error_vs_true_frequency_bin, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_error_vs_true_frequency_bin.png")
```

![](../figures/neutral_error_vs_true_frequency_bin.png)<!-- -->

## Check the SNPs with highest error

``` r
filter(joined_frequency_final, coverage==8, sample_size==160) %>%
  mutate(error=frequency-estimated_frequency) %>%
  arrange(desc(abs(error))) %>%
  head(n=20)
```

    ## # A tibble: 20 x 13
    ##    type  position frequency base  major minor anc   estimated_frequ…  nInd
    ##    <chr>    <dbl>     <dbl> <chr> <chr> <chr> <chr>            <dbl> <dbl>
    ##  1 m1    23704228     0.683 C     G     C     G                0.567   160
    ##  2 m1    28110618     0.36  T     G     T     G                0.470   160
    ##  3 m1    28604410     0.469 T     G     T     G                0.366   160
    ##  4 m1    19794246     0.567 A     G     A     G                0.669   160
    ##  5 m1    19462710     0.496 G     T     G     T                0.395   160
    ##  6 m1    23705057     0.464 T     G     T     G                0.364   160
    ##  7 m1    23705146     0.458 C     T     C     T                0.358   160
    ##  8 m1    19794252     0.414 C     T     C     T                0.316   160
    ##  9 m1           6     0.674 T     C     T     C                0.773    22
    ## 10 m1    28110381     0.359 C     T     C     T                0.458   160
    ## 11 m1    19794182     0.584 C     A     C     A                0.681   160
    ## 12 m1    23543799     0.352 C     A     C     A                0.256   160
    ## 13 m1      824660     0.337 C     T     C     T                0.433   160
    ## 14 m1    27387986     0.487 T     A     T     A                0.583   160
    ## 15 m1    18417553     0.523 A     C     A     C                0.427   160
    ## 16 m1    11118911     0.527 G     T     G     T                0.431   160
    ## 17 m1    24595214     0.376 G     C     G     C                0.280   160
    ## 18 m1    24354179     0.494 T     C     T     C                0.589   160
    ## 19 m1    24595222     0.379 G     T     G     T                0.285   160
    ## 20 m1    23775604     0.387 T     G     T     G                0.481   160
    ## # … with 4 more variables: coverage <dbl>, sample_size <dbl>,
    ## #   frequency_bin <fct>, error <dbl>

## True frequency distribution of false negatives

``` r
false_negatives_final_count <- count(false_negatives_final, coverage, sample_size)
ggplot(false_negatives_final, aes(x=frequency)) +
  geom_histogram() +
  geom_text(data=false_negatives_final_count, x=0.3, y=200000, aes(label=paste0("n=", n))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Esimated frequency distribution of false positives

``` r
false_positives_final_count <- count(false_positives_final, coverage, sample_size)
ggplot(false_positives_final, aes(x=estimated_frequency)) +
  geom_histogram() +
  geom_text(data=false_positives_final_count, x=0.8, y=1500, aes(label=paste0("n=", n))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

## Read windowed thetas estimated by `realSFS`

``` r
i=1
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40, 80, 160)){
    ## read in estimated maf
    thetas <- read_tsv(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.thetas.idx.pestPG")) %>%
      transmute(position=WinCenter, theta_w=tW/nSites, theta_t=tP/nSites, tajima_d=Tajima) %>%
      mutate(coverage=coverage, sample_size=sample_size) %>%
      gather(key=summary_stats, value=value, 2:4)
    average_thetas <- read_tsv(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.average_thetas.idx.pestPG")) %>%
      transmute(theta_w=tW/nSites, theta_t=tP/nSites, tajima_d=Tajima) %>%
      mutate(coverage=coverage, sample_size=sample_size)
    ## compile the final files for plotting
    if (i==1){
      thetas_final <- thetas
      average_thetas_final <- average_thetas
    } else {
      thetas_final <- bind_rows(thetas_final, thetas)
      average_thetas_final <- bind_rows(average_thetas_final, average_thetas)
    }
    i=i+1
  }
}
```

## real theta values with the entire population

``` r
real_theta_t <- sum(2*mutations_final$frequency*(1-mutations_final$frequency))/30000000
real_theta_w <- dim(mutations_final)[1]/(30000000*sum(1/(1:999)))
real_tajima_d <- real_theta_t-real_theta_w
tibble(real_theta_t=round(real_theta_t, 5), real_theta_w=round(real_theta_w, 5)) %>%
  kable()
```

| real\_theta\_t | real\_theta\_w |
| -------------: | -------------: |
|        0.00394 |        0.00433 |

## thetas estimated from the estimated allele frequencies (not from `realSFS`)

Tajima’s estimator

``` r
group_by(joined_frequency_final, coverage, sample_size) %>%
  summarise(theta_t=round(sum(2*estimated_frequency*(1-estimated_frequency))/30000000,5)) %>%
  ungroup() %>%
  spread(key = sample_size, value = theta_t) %>%
  kable()
```

| coverage |       5 |      10 |      20 |      40 |      80 |     160 |
| -------: | ------: | ------: | ------: | ------: | ------: | ------: |
|     0.25 | 0.00013 | 0.00055 | 0.00149 | 0.00245 | 0.00301 | 0.00343 |
|     0.50 | 0.00055 | 0.00148 | 0.00243 | 0.00301 | 0.00343 | 0.00366 |
|     1.00 | 0.00146 | 0.00242 | 0.00301 | 0.00342 | 0.00365 | 0.00377 |
|     2.00 | 0.00239 | 0.00301 | 0.00341 | 0.00365 | 0.00376 | 0.00382 |
|     4.00 | 0.00299 | 0.00341 | 0.00364 | 0.00377 | 0.00382 | 0.00385 |
|     8.00 | 0.00335 | 0.00363 | 0.00376 | 0.00383 | 0.00386 | 0.00387 |

Watterson’s estimator

``` r
group_by(joined_frequency_final, coverage, sample_size) %>%
  summarise(theta_w=round(n()/30000000/sum(1/(1:(unique(sample_size)-1))), 5)) %>%
  spread(key = sample_size, value = theta_w) %>%
  kable()
```

| coverage |       5 |      10 |      20 |      40 |      80 |     160 |
| -------: | ------: | ------: | ------: | ------: | ------: | ------: |
|     0.25 | 0.00014 | 0.00044 | 0.00098 | 0.00143 | 0.00166 | 0.00190 |
|     0.50 | 0.00058 | 0.00122 | 0.00171 | 0.00194 | 0.00218 | 0.00235 |
|     1.00 | 0.00162 | 0.00213 | 0.00235 | 0.00254 | 0.00269 | 0.00277 |
|     2.00 | 0.00284 | 0.00296 | 0.00306 | 0.00315 | 0.00320 | 0.00323 |
|     4.00 | 0.00393 | 0.00383 | 0.00378 | 0.00376 | 0.00374 | 0.00372 |
|     8.00 | 0.00483 | 0.00454 | 0.00439 | 0.00429 | 0.00422 | 0.00415 |

## thetas estimated from `realSFS`

Tajima’s estimator

``` r
select(average_thetas_final, theta_t, coverage, sample_size) %>%
  spread(key = sample_size, value=theta_t) %>%
  kable()
```

| coverage |         5 |        10 |        20 |        40 |        80 |       160 |
| -------: | --------: | --------: | --------: | --------: | --------: | --------: |
|     0.25 | 0.0034168 | 0.0031124 | 0.0032314 | 0.0034008 | 0.0035915 | 0.0037341 |
|     0.50 | 0.0031246 | 0.0032566 | 0.0034224 | 0.0036059 | 0.0037399 | 0.0038179 |
|     1.00 | 0.0032922 | 0.0034674 | 0.0036276 | 0.0037527 | 0.0038253 | 0.0038660 |
|     2.00 | 0.0035273 | 0.0036772 | 0.0037772 | 0.0038415 | 0.0038762 | 0.0038965 |
|     4.00 | 0.0037369 | 0.0038194 | 0.0038670 | 0.0038948 | 0.0039081 | 0.0039118 |
|     8.00 | 0.0038659 | 0.0039005 | 0.0039167 | 0.0039249 | 0.0039249 | 0.0039203 |

Watterson’s estimator

``` r
select(average_thetas_final, theta_w, coverage, sample_size) %>%
  spread(key = sample_size, value=theta_w) %>%
  kable()
```

| coverage |         5 |        10 |        20 |        40 |        80 |       160 |
| -------: | --------: | --------: | --------: | --------: | --------: | --------: |
|     0.25 | 0.0042660 | 0.0020508 | 0.0019516 | 0.0018931 | 0.0019744 | 0.0020936 |
|     0.50 | 0.0024890 | 0.0022838 | 0.0022244 | 0.0022713 | 0.0023722 | 0.0024618 |
|     1.00 | 0.0027675 | 0.0026765 | 0.0026811 | 0.0027483 | 0.0028212 | 0.0028751 |
|     2.00 | 0.0032025 | 0.0032222 | 0.0032658 | 0.0033009 | 0.0033326 | 0.0033450 |
|     4.00 | 0.0035851 | 0.0036249 | 0.0036617 | 0.0036731 | 0.0036847 | 0.0036793 |
|     8.00 | 0.0038146 | 0.0038444 | 0.0038651 | 0.0038656 | 0.0038690 | 0.0038578 |

## Plot Watterson’s estimator and Tajima’s estimator of theta and Tajima’s D in 10,000bp fixed windows

``` r
filter(thetas_final, summary_stats !="tajima_d") %>%
  ggplot(aes(x=position, y=value, color=summary_stats)) +
  geom_line(size=0.2) +
  # geom_text(data=false_positives_final_count, x=0.8, y=1500, aes(label=paste0("n=", n))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
filter(thetas_final, summary_stats =="tajima_d") %>%
  ggplot(aes(x=position, y=value, color=summary_stats)) +
  geom_line(size=0.2) +
  geom_hline(yintercept = 0, color="black") +
  # geom_text(data=false_positives_final_count, x=0.8, y=1500, aes(label=paste0("n=", n))) +
  facet_grid(coverage~sample_size) +
  theme_cowplot()
```

![](data_analysis_neutral_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

I will annonate each figure with the chromosome average statistics later
on.

# Compare individual barcoding with Pool-seq

I assumed equal sequencing output from all individuals in this step.
This is highly unlikely when no individual barcode is given, but it
represents the best case scenario for Pool-seq.

A better way to do this is to merge all the 20x bam files from all
individuals, and subsample from this huge merged bam file. This can be
costly computationally, but can still be feasible if we want to go with
it.

I used the same set of SNPs obtained from ANGSD and calculated allele
frequency based on total allele count across the population. This is to
make the results from the two methods more comparable. In reality,
however, the SNP calling result will be different if no individual
barcode is provided.

## Get minor allele frequencies estimated from Pool-seq

``` r
i=1
for (coverage in c(0.25,0.5,1,2,4,8)){
  for (sample_size in c(5,10,20,40, 80, 160)){
    ## read in allele count
    allele_count <- read_tsv(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.pos.gz")) %>%
      bind_cols(read_tsv(paste0("../neutral_sim/rep_1/angsd/bam_list_", sample_size, "_", coverage, "x.counts.gz"))) %>%
      mutate(ancestral_allele = ancestral$ancestral[pos]) %>%
    ## turn allele count into maf
      rowwise() %>%
      mutate(maf = count_to_maf(ancestral_allele, totA, totC, totG, totT), position=pos) %>%
      ungroup()
    ## join estimated maf with true snps and only keep the snps that exist in estimated maf
    joined_frequency <- right_join(mutations_final, allele_count, by="position") %>%
      select(-ancestral) %>%
      mutate(coverage=coverage, sample_size=sample_size, 
             frequency=ifelse(is.na(frequency) & !(position %in% substitutions_pos), 0, frequency),
             frequency=ifelse(is.na(frequency) & position %in% substitutions_pos, 1, frequency))
    ## compile the final files for plotting
    if (i==1){
      joined_frequency_final <- joined_frequency
    } else {
      joined_frequency_final <- bind_rows(joined_frequency_final, joined_frequency)
    }
    i=i+1
  }
}
write_tsv(joined_frequency_final, "../neutral_sim/rep_1/angsd/joined_frequency_final_poolseq.tsv")
```

``` r
joined_frequency_final <- read_tsv("../neutral_sim/rep_1/angsd/joined_frequency_final_poolseq.tsv") %>%
  mutate(estimated_frequency=maf, frequency_bin = cut(frequency, breaks = 0:10/10, include.lowest = T), error=estimated_frequency-frequency) 
joined_summary <- summarise_by_design(joined_frequency_final)
```

## Plot estimated allele frequency vs. true allele frequency (this includes the false positives but not the false negatives)

``` r
estimated_vs_true_frequency_pool <- plot_frequency(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_estimated_vs_true_frequency_pool.png", estimated_vs_true_frequency_pool, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_estimated_vs_true_frequency_pool.png")
```

![](../figures/neutral_estimated_vs_true_frequency_pool.png)<!-- -->

## Plot estimated allele frequency vs. true allele frequency in bins (this includes the false positives but not the false negatives)

``` r
estimated_vs_true_frequency_pool_bin <- plot_frequency_in_bins(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_estimated_vs_true_frequency_pool_bin.png", estimated_vs_true_frequency_pool_bin, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_estimated_vs_true_frequency_pool_bin.png")
```

![](../figures/neutral_estimated_vs_true_frequency_pool_bin.png)<!-- -->

## Plot absolute values of error vs. true allele frequency in bins (this includes the false positives but not the false negatives)

``` r
error_vs_true_frequency_pool_bin <- plot_error_in_bins(joined_frequency_final, joined_summary)
ggsave("../figures/neutral_error_vs_true_frequency_pool_bin.png", error_vs_true_frequency_pool_bin, height = 12, width=15, units = "in")
```

``` r
include_graphics("../figures/neutral_error_vs_true_frequency_pool_bin.png")
```

![](../figures/neutral_error_vs_true_frequency_pool_bin.png)<!-- -->
