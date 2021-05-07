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
get_neutrality_stats <- function(x, subdir="angsd_gatk"){
  i=1
  for (coverage in c(0.25,0.5,1,2,4,8)){
    for (sample_size in c(5,10,20,40,80, 160)){
      ## read in estimated fst
      neutrality_stats <- read_tsv(paste0(x, subdir, "/bam_list_p1_", sample_size, "_", coverage, "x_all_sites.windowed_thetas.idx.pestPG")) %>%
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