library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)
library(broom)
library(GGally)

#Written for the analysis of integrated library expression at 8 µM. 2 biological
#replicates and 2 technical replicates each for DNA amplification
#D##_BC: DNA biological replicate # technical replicate #
#R#8_BC: RNA biological replicate # at 8 µM


#Load index and bcmap files-----------------------------------------------------

bc_DNA_1_1 <- read_tsv('BCreads_txts/D11_BC.txt')
bc_DNA_1_2 <- read_tsv('BCreads_txts/D12_BC.txt')
bc_DNA_2_1 <- read_tsv('BCreads_txts/D21_BC.txt')
bc_DNA_2_2 <- read_tsv('BCreads_txts/D22_BC.txt')
bc_RNA_1 <- read_tsv('BCreads_txts/R18_BC.txt')
bc_RNA_2 <- read_tsv('BCreads_txts/R28_BC.txt')


#Load barcode mapping table, remember sequences are rcomp due to sequencing 
#format

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff)


#pick out SP3 and SP5 in the bcmap that were used in this assay

SP3_SP5_map <- barcode_map %>%
  mutate(subpool = ifelse(
    startsWith(name, 'subpool'), 
    substr(name, 1, 8), 
    'control')) %>%
  filter(subpool != 'subpool2') %>%
  filter(subpool != 'subpool4')


#Join reads to bcmap------------------------------------------------------------

#Join BC reads to BC mapping, keeping the reads only appearing in barcode 
#mapping and replacing na with 0 reads.

bc_map_join_bc <- function(df1, df2) {
  df2 <- df2 %>%
    mutate(norm = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(norm = if_else(is.na(norm), 
                                0, 
                                norm)) %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                               as.integer(0), 
                               num_reads))
  return(keep_bc)
}

bc_join_DNA_1_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_1)
bc_join_DNA_1_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_2)
bc_join_DNA_2_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_1)
bc_join_DNA_2_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_2)
bc_join_RNA_1 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_1)
bc_join_RNA_2 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_2)


#Med analysis-------------------------------------------------------------------

#Join DNA BC reads > 6, take average norm. reads and join to RNA. Take ratio of 
#RNA/DNA norm reads

ave_dna_join_rna_rep <- function(df1, df2, df3) {
  filter_reads_1 <- filter(df1, num_reads > 6)
  filter_reads_2 <- filter(df2, num_reads > 6)
  DNA_join <- inner_join(filter_reads_1, filter_reads_2, 
                         by = c("barcode", "name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")) %>%
    mutate(ave_norm_DNA = (norm_DNA_tr1 + norm_DNA_tr2)/2)
  DNA_RNA_join <- left_join(DNA_join, df3,
                             by = c("barcode", "name", "subpool", 
                                    "most_common")) %>%
    rename(num_reads_RNA = num_reads) %>%
    rename(norm_RNA = norm) %>%
    mutate(ratio = norm_RNA/ave_norm_DNA)
  print('processed dfs in order of (DNA tr1, DNA tr2, RNA) in 
        bc_dna_biol_rep(df1, df2, df3)')
  return(DNA_RNA_join)
}

bc_ave_DNA_RNA_1 <- ave_dna_join_rna_rep(bc_join_DNA_1_1, bc_join_DNA_1_2, 
                                         bc_join_RNA_1)
bc_ave_DNA_RNA_2 <- ave_dna_join_rna_rep(bc_join_DNA_2_1, bc_join_DNA_2_2, 
                                         bc_join_RNA_2)


#Count barcodes per variant per DNA and RNA, set minimum of 8 BC's per variant 
#in both samples, take median RNA/DNA per variant, then per variant determine 
#the median absolute deviation of all barcode ratios. Then filter out variants 
#with 0 median expression

ratio_bc_med_var <- function(df) {
  bc_count_DNA <- df %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n()) %>%
    filter(barcodes_RNA > 7)
  bc_DNA_RNA <- inner_join(bc_count_DNA, bc_count_RNA, 
                           by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  bc_min_8_df <- left_join(bc_DNA_RNA, df, 
                           by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  med_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio))
  mad_ratio <- bc_min_8_df %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = mad(ratio))
  med_mad <- inner_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio)) %>%
    mutate(mad_over_med = if_else(
      is.na(mad_over_med),
      as.double(0), 
      mad_over_med))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup() %>%
    filter(med_ratio > 0)
  return(bc_med)
}

med_ratio_1 <- ratio_bc_med_var(bc_ave_DNA_RNA_1)
med_ratio_2 <- ratio_bc_med_var(bc_ave_DNA_RNA_2)


#combine biological replicates

rep_1_2 <- inner_join(med_ratio_1, med_ratio_2,
             by = c("name", "subpool", "most_common"),
             suffix = c('_br1', '_br2'))

output_int <- rep_1_2 %>%
  write.table(
    "rep_1_2.txt", 
    sep = '\t', row.names = FALSE)

#plot replicability, output table and correlation

rep_1_2_log10 <- var_log10(rep_1_2)

p_var_med_rep <- rep_1_2 %>%
  ggplot(aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.5) + 
  geom_point(data = filter(rep_1_2, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.001, 20), breaks = c(0.01, 0.1, 1, 10)) + 
  scale_y_log10(limits = c(0.001, 20), breaks = c(0.01, 0.1, 1, 10)) +
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) 

save_plot('plots/p_var_med_rep.pdf', p_var_med_rep, scale = 1.3,
          base_width = 2.7, base_height = 2.35)

pearsons_med <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(rep_1_2_log10$med_ratio_br1, 
                         rep_1_2_log10$med_ratio_br2, 
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3),
               round(cor(filter(rep_1_2_log10, 
                                subpool == 'subpool3')$med_ratio_br1,
                         filter(rep_1_2_log10, 
                                subpool == 'subpool3')$med_ratio_br2,
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3),
               round(cor(filter(rep_1_2_log10, 
                                subpool == 'subpool5')$med_ratio_br1,
                         filter(rep_1_2_log10, 
                                subpool == 'subpool5')$med_ratio_br2,
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3)))

write_csv(pearsons_med, 'pearsons_med.csv')


#After combining, rename backgrounds to simplified names, make background column 
#(excluding controls), separate out background values in each dataset and left 
#join to original dataset. Normalize expression of each variant to its 
#background in that biological replicate. Determine average normalized 
#expression across biological replicates.

back_norm <- function(df1) {
  gsub_1_2 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, nchar(background)-5, nchar(background)))
  backgrounds <- gsub_1_2 %>%
    filter(startsWith(name, 'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, med_ratio_br1, med_ratio_br2) %>%
    rename(med_ratio_br1_back = med_ratio_br1) %>%
    rename(med_ratio_br2_back = med_ratio_br2) 
  back_join_norm <- left_join(gsub_1_2, backgrounds, by = 'background') %>%
    mutate(med_ratio_br1_norm = med_ratio_br1/med_ratio_br1_back) %>%
    mutate(med_ratio_br2_norm = med_ratio_br2/med_ratio_br2_back) %>%
    mutate(ave_med_ratio_norm = (med_ratio_br1_norm + med_ratio_br2_norm)/2)
}

int_back_norm_rep_1_2 <- back_norm(rep_1_2)

int_back_norm_pc_spGl4 <- rep_1_2 %>%
  filter(name == 'pGL4.29 Promega 1-63 + 1-87' | name == 'pGL4.29 Promega 1-87') %>%
  mutate(name = str_c(name, '_scramble pGL4.29 Promega 1-63 + 1-87')) %>%
  mutate(subpool = 'subpool3') %>%
  rbind(rep_1_2) %>%
  back_norm()


#determine the log(RNA/DNA) for each sample

var_log2 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, 
              funs(log2(.))
    )
  return(log_ratio_df)
}

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, 
              funs(log10(.))
    )
  return(log_ratio_df)
}

rep_1_2_log10 <- var_log10(rep_1_2)

int_back_norm_rep_1_2_log10 <- var_log10(int_back_norm_rep_1_2)

int_back_norm_pc_spGl4_log10 <- var_log10(int_back_norm_pc_spGl4)


#Sum analysis-------------------------------------------------------------------

#sum unique barcodes and normalized bc reads across barcodes per variant. Set-up
#a minimum of 8 barcodes per sample to reliably sum. 

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n()) %>%
    filter(barcodes > 7)
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = norm) %>%
    rename(sum = n)
  bc_sum <- right_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA_1_1 <- var_sum_bc_num(bc_join_DNA_1_1)
variant_counts_DNA_1_2 <- var_sum_bc_num(bc_join_DNA_1_2)
variant_counts_DNA_2_1 <- var_sum_bc_num(bc_join_DNA_2_1)
variant_counts_DNA_2_2 <- var_sum_bc_num(bc_join_DNA_2_2)
variant_counts_RNA_1 <- var_sum_bc_num(bc_join_RNA_1)
variant_counts_RNA_2 <- var_sum_bc_num(bc_join_RNA_2)


#Average DNA variant sum across technical replicates. Left_join with RNA that is
#the same biological replicate and determine variant expression as = sum nrpm 
#RNA/sum nrpm DNA. Filter out expression = 0

ave_dna_join_rna_var <- function(df1, df2, df3) {
  DNA_join <- inner_join(df1, df2, 
                         by = c("name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")) %>%
    mutate(ave_sum_DNA = (sum_DNA_tr1 + sum_DNA_tr2)/2)
  DNA_RNA_join <- left_join(DNA_join, df3,
                            by = c("name", "subpool", "most_common")) %>%
    rename(sum_RNA = sum) %>% 
    rename(barcodes_RNA = barcodes) %>%
    filter(sum_RNA > 0) %>%
    mutate(ratio = sum_RNA/ave_sum_DNA)
  print('processed dfs in order of (DNA tr1, DNA tr2, RNA) in 
        ave_dna_join_rna_var(df1, df2, df3)')
  return(DNA_RNA_join)
}

RNA_DNA_br1 <- ave_dna_join_rna_var(variant_counts_DNA_1_1, 
                                    variant_counts_DNA_1_2, 
                                    variant_counts_RNA_1)
RNA_DNA_br2 <- ave_dna_join_rna_var(variant_counts_DNA_2_1, 
                                    variant_counts_DNA_2_2, 
                                    variant_counts_RNA_2)

rep_1_2_sum <- inner_join(RNA_DNA_br1, RNA_DNA_br2, 
                          by = c('name', 'subpool', 'most_common'),
                          suffix = c('_br1', '_br2'))


#Write output table

output_int <- rep_1_2_sum %>%
  write.table(
    "rep_1_2_sum.txt", 
    sep = '\t', row.names = FALSE)


#Determine log expression

rep_1_2_sum_log10 <- var_log10(rep_1_2_sum)


#Sum replicability

p_var_sum_rep <- rep_1_2_sum %>%
  ggplot(aes(ratio_br1, ratio_br2)) +
  geom_point(alpha = 0.1, size = 0.75) +
  geom_point(data = filter(rep_1_2_sum, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.5) + 
  geom_point(data = filter(rep_1_2_sum, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.004, 20)) + 
  scale_y_log10(limits = c(0.004, 20)) +
  background_grid(major = 'xy', minor = 'none') + 
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) 

save_plot('plots/p_var_sum_rep.pdf', p_var_sum_rep, scale = 1.3,
          base_width = 2.5, base_height = 2.35)

pearsons_sample <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(rep_1_2_sum_log10$ratio_br1, 
                         rep_1_2_sum_log10$ratio_br2, 
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3),
               round(cor(filter(rep_1_2_sum_log10, 
                                subpool == 'subpool3')$ratio_br1,
                         filter(rep_1_2_sum_log10, 
                                subpool == 'subpool3')$ratio_br2,
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3),
               round(cor(filter(rep_1_2_sum_log10, 
                                subpool == 'subpool5')$ratio_br1,
                         filter(rep_1_2_sum_log10, 
                                subpool == 'subpool5')$ratio_br2,
                         use = "pairwise.complete.obs", 
                         method = "pearson"), 3)))

write_csv(pearsons_sample, 'pearsons_sample.csv')


#Separate into subpools----------------------------------------------------------------------

#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) that vary in 
#distance from one another by 0 (no inner flanks), 5, 10, 15, 20 and 70 bp (all but 0 appear 
#as -4 bp spacing). Each site distance combination is then moved along the backgrounds at 1 
#bp increments starting from closest to the minP. Separation lists the spacing between sites 
#and distance (start of consensus and flanks). Added 2 to all distances to measure to start of
#BS and not to flank. Added 4 to all spacing but 0 to measure difference between start of 
#sites. Also took average of log2 med BC expression between biological replicates for plotting

subpool3 <- 
  filter(int_back_norm_rep_1_2_log10, subpool == "subpool3") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, 
           into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "fluff4"),
           sep = "_", convert = TRUE
           ) %>%
  select(-subpool, -fluff2, -fluff3, -fluff4) %>%
  mutate(dist = as.integer(dist + 2)) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  as.integer(spacing + 4), as.integer(spacing))) 

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from furthest to 
#the minP. These sites are filled with sites of either the consensus site, a weak site or no 
#site. Both the weak and consensus sites are flanked by the same flanking sequence. 

subpool5 <- 
  filter(int_back_norm_rep_1_2, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", 
    "fluff"), sep = "_") %>%
  select(-subpool, -fluff) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus")) %>%
  mutate(weak = str_detect(site1, "weak") +
           str_detect(site2, "weak") +
           str_detect(site3, "weak") +
           str_detect(site4, "weak") +
           str_detect(site5, "weak") +
           str_detect(site6, "weak")) %>%
  mutate(nosite = str_detect(site1, "nosite") +
           str_detect(site2, "nosite") +
           str_detect(site3, "nosite") +
           str_detect(site4, "nosite") +
           str_detect(site5, "nosite") +
           str_detect(site6, "nosite")) %>%
  mutate(total_sites = consensus + weak) %>%
  mutate(site_combo = 
           ifelse(weak == 0 & consensus > 0, 
                  'consensus', 'mixed')) %>%
  mutate(site_type = 
           ifelse(consensus == 0 & weak > 0, 
                  'weak', site_combo))

controls <- 
  filter(rep_1_2, subpool == "control") %>%
  ungroup() %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)


#Write table to compare expression to other analyses-----------------------------------------

output_int <- rep_1_2 %>%
  write.table(
    "rep_1_2.txt", 
    sep = '\t', row.names = FALSE)

output_int_controls <- controls %>%
  write.table(
    'int_controls_rep_1_2.txt',
    sep = '\t', row.names = FALSE
  )


#plot variant statistics and rep plots-------------------------------------------------------

#replicates

p_var_med_rep <- rep_1_2 %>%
  ggplot(aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(alpha = 0.1, size = 1) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             fill = 'orange', shape = 21, size = 1.75) + 
  geom_point(data = filter(rep_1_2, name == 'pGL4.29 Promega 1-63 + 1-87'), 
             fill = 'red', shape = 21, size = 1.75) +
  geom_abline(slope = 1, color = 'grey60') +
  annotation_logticks(scaled = TRUE) +
  xlab("Expression (a.u.) replicate 1") +
  ylab("Expression (a.u.) replicate 2") +
  scale_x_log10(limits = c(0.001, 15)) + 
  scale_y_log10(limits = c(0.001, 15)) +
  background_grid(major = 'xy', minor = 'none') + 
  theme(strip.background = element_rect(colour="black", fill="white"),
        axis.line.y = element_line(), panel.spacing.x=unit(1, "lines")) 

save_plot('plots/p_var_med_rep_line_sp5.pdf', p_var_med_rep, scale = 1.3,
          base_width = 2.5, base_height = 2.35)

pearsons_sample <- tibble(
  sample = c('all', 'subpool3', 'subpool5'),
  pearsons = c(round(cor(rep_1_2_log10$med_ratio_br1, 
                         rep_1_2_log10$med_ratio_br2, 
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(rep_1_2_log10, 
                                subpool == 'subpool3')$med_ratio_br1,
                         filter(rep_1_2_log10, 
                                subpool == 'subpool3')$med_ratio_br2,
                         use = "pairwise.complete.obs", method = "pearson"), 3),
               round(cor(filter(rep_1_2_log10, 
                                subpool == 'subpool5')$med_ratio_br1,
                         filter(rep_1_2_log10, 
                                subpool == 'subpool5')$med_ratio_br2,
                         use = "pairwise.complete.obs", method = "pearson"), 3)))

write_csv(pearsons_sample, 'pearsons_sample.csv')


#Mad over median analysis

med_mad_gather <- function(df1) {
  med <- df1 %>%
    select(subpool, name, med_ratio_br1, med_ratio_br2) %>%
    gather(med_ratio_br1, med_ratio_br2, key = condition, value = med) %>%
    separate(condition, into = c("fluff1", "fluff2", "condition"),
             sep = "_", convert = TRUE) %>% 
    select(-fluff1, -fluff2)
  mad_over_med <- df1 %>%
    select(subpool, name, mad_over_med_br1, mad_over_med_br2) %>%
    gather(mad_over_med_br1, mad_over_med_br2, 
           key = condition, value = mad_over_med) %>%
    separate(condition, into = c("fluff1", "fluff2", "fluff3", "condition"),
             sep = "_", convert = TRUE) %>% 
    select(-fluff1, -fluff2, -fluff3)
  med_mad <- inner_join(med, mad_over_med, 
                        by = c('subpool', 'name', 'condition'))
  return(med_mad)
}

med_mad_rep_1_2 <- med_mad_gather(rep_1_2)

ggplot(med_mad_rep_1_2, aes(mad_over_med)) +
  geom_density(kernel = 'gaussian') +
  facet_grid(subpool ~ condition) + 
  panel_border() +
  xlab('MAD/Med')

ggplot(med_mad_rep_1_2, aes(med, mad_over_med)) +
  ylab('Mad/Med') +
  xlab('Median') +
  geom_point(alpha = 0.3)

rep_1_2_log10_madomed <- rep_1_2_0rm(rep_1_2) %>%
  filter(mad_over_med_br1 == 1.48260000 | mad_over_med_br2 == 1.48260000) %>%
  var_log10()

p_var_med_ratio_madmed1 <- ggplot(NULL, aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(data = rep_1_2_log10, color = 'black', alpha = 0.2) + 
  geom_point(data = rep_1_2_log10_madomed, color = 'red', alpha = 0.2) + 
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant median\n RNA/DNA BR 1") +
  ylab("log10 variant median\n RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) + 
  scale_y_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) +
  annotate("text", x = -2, y = 0.5,
           label = paste('rtot =', 
                         round(cor(rep_1_2_log10$med_ratio_br1, 
                                   rep_1_2_log10$med_ratio_br2,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

save_plot('plots/p_var_med_ratio_madmed1.png', p_var_med_ratio_madmed1)


#Plot subpool expression features-----------------------------------------------

#Subpool 3

p_subpool3_spa_back_norm <- ggplot(
  filter(subpool3, background != 'v chr9'), 
  aes(x = dist)
  ) + 
  geom_point(aes(y = med_ratio_br1_norm), alpha = 0.7, color = '#287D8EFF') +
  geom_point(aes(y = med_ratio_br2_norm), alpha = 0.7, color = '#95D840FF') +
  geom_smooth(aes(y = ave_med_ratio_norm), color = '#482677FF', span = 0.13) +
  facet_grid(spacing ~ background) + 
  ylab('log10 average normalized median BC expression') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
    breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/p_subpool3_spa_back_norm.png', p_subpool3_spa_back_norm, 
          base_width = 46, base_height = 17, scale = 0.35)

#Subpool 5

p_subpool5_cons_mix_weak_sitenum <- ggplot(
  NULL, aes(as.factor(total_sites), ave_med_ratio_norm)
  ) + 
  facet_wrap(background ~ site_type) + 
  geom_boxplot(data = filter(subpool5, total_sites > 0 & total_sites < 6)) +
  geom_boxplot(data = filter(subpool5, site_type == 'mixed' & total_sites == 6)) +
  geom_point(data = filter(subpool5, site_type != 'mixed' & total_sites == 6)) +
  theme(panel.spacing.x=unit(0, "lines")) +
  annotation_logticks(sides = 'l') +
  background_grid(major = 'y') +
  panel_border() + 
  xlab("Total Number of Binding Sites") +
  scale_y_continuous("log10 median BC expression")

p_subpool5_cons_mix_sitenum <- ggplot(NULL, aes(site_type, ave_med_ratio_norm)) +
  facet_grid(background ~ total_sites) + 
  geom_boxplot(data = filter(subpool5, total_sites > 0 & total_sites < 6)) +
  geom_boxplot(data = filter(subpool5, site_type == 'mixed' & total_sites == 6)) +
  geom_point(data = filter(subpool5, site_type != 'mixed' & total_sites == 6)) +
  theme(panel.spacing.x=unit(0, "lines"), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  background_grid(major = 'y') +
  panel_border() + 
  xlab("Binding site type") +
  scale_y_continuous("log10 median BC expression") +
  annotation_logticks(sides = 'l')

save_plot('plots/p_subpool5_cons_mix_sitenum.png', p_subpool5_cons_mix_sitenum, 
          base_width = 46, base_height = 17, scale = 0.35)


#Plot BC statistics and rep plots------------------------------------------------------------

#BC num reads + #BC's

p_BC_num_reads_viol_full <- ggplot(NULL) +
  geom_jitter(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR1', y = num_reads_DNA_tr1), size = 0.5, alpha = 0.1) +
  geom_jitter(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR2', y = num_reads_DNA_tr2), size = 0.5, alpha = 0.1) +
  geom_jitter(data = bc_ave_DNA_RNA_1, 
              aes(x = 'RNA BR1', y = num_reads_RNA), size = 0.5, alpha = 0.1) +
  geom_jitter(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR1', y = num_reads_DNA_tr1), size = 0.5, alpha = 0.1) +
  geom_jitter(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR2', y = num_reads_DNA_tr2), size = 0.5, alpha = 0.1) +
  geom_jitter(data = bc_ave_DNA_RNA_2, 
              aes(x = 'RNA BR2', y = num_reads_RNA), size = 0.5, alpha = 0.1) +
  xlab("") +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  ylab("Reads per BC")
  
save_plot('plots/BC_num_reads_viol_full.png', p_BC_num_reads_viol_full, 
          scale = 1.3, base_width = 6.5, base_height = 3)

p_BC_num_reads_box_zoom <- ggplot(NULL) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR1', y = num_reads_DNA_tr1)) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR2', y = num_reads_DNA_tr2)) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'RNA BR1', y = num_reads_RNA)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR1', y = num_reads_DNA_tr1)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR2', y = num_reads_DNA_tr2)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'RNA BR2', y = num_reads_RNA)) +
  xlab("") +
  ylab("Reads per BC") + ylim(0, 25)

save_plot('plots/BC_num_reads_box_zoom.png', p_BC_num_reads_box_zoom, 
          scale = 1.3,base_width = 6.5, base_height = 3)

p_BC_number <- ggplot(NULL) +
  geom_boxplot(data = rep_1_2, aes(x = 'DNA_br1', y = barcodes_DNA_br1)) +
  geom_boxplot(data = rep_1_2, aes(x = 'DNA_br2', y = barcodes_DNA_br2)) +
  geom_boxplot(data = rep_1_2, aes(x = 'RNA_br1', y = barcodes_RNA_br1)) +
  geom_boxplot(data = rep_1_2, aes(x = 'RNA_br2', y = barcodes_RNA_br2)) +
  xlab('') +
  ylab('BCs/variant') +
  background_grid(major = 'y')

save_plot('plots/BC_number.png', p_BC_number, scale = 1.3)


#BC normalized reads rep (need to log transform normalized reads first to compare, RNA and 
#ave DNA comparison required merging BC's that are present in each br (DNA and RNA))

p_bc_rep_DNA1 <- ggplot(data = var_log10(bc_ave_DNA_RNA_1), 
                        aes(norm_DNA_tr1, norm_DNA_tr2)) +
  geom_point(alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads TR 1") +
  ylab("Log10 norm. BC reads TR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5,2.5)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5,2.5)) + 
  annotate("text", x = 0, y = 2, 
           label = paste('r =', 
                         round(cor(var_log10(bc_ave_DNA_RNA_1)$norm_DNA_tr1, 
                                   var_log10(bc_ave_DNA_RNA_1)$norm_DNA_tr2,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_bc_rep_DNA2 <- ggplot(data = var_log10(bc_ave_DNA_RNA_2), 
                        aes(norm_DNA_tr1, norm_DNA_tr2)) +
  geom_point(alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads TR 1") +
  ylab("Log10 norm. BC reads TR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5,2.5)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5,2.5)) + 
  annotate("text", x = 0, y = 2, 
           label = paste('r =', 
                         round(cor(var_log10(bc_ave_DNA_RNA_2)$norm_DNA_tr1, 
                                   var_log10(bc_ave_DNA_RNA_2)$norm_DNA_tr2,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_bc_rep_grid <- plot_grid(p_bc_rep_DNA1, p_bc_rep_DNA2, 
                           labels = c('BR1 DNA', 'BR2 DNA'), align = 'h', 
                           hjust = -1.2, vjust = 0.1, scale = 0.9)

save_plot('plots/p_bc_rep_grid.png', p_bc_rep_grid, 
          base_height = 3, base_width = 6, scale = 1.3)


#Fitting a model to the 6-site library------------------------------------------

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}


#Linear models

ind_site_ind_back <- function(df) {
  model <- lm(ave_med_ratio_norm ~ site1 + site2 + site3 + site4 + site5 + site6 + background, 
              data = df)
}

subpool5_ncw <- subpool5 %>%
  mutate(site1 = gsub('nosite', 'anosite', site1)) %>%
  mutate(site2 = gsub('nosite', 'anosite', site2)) %>%
  mutate(site3 = gsub('nosite', 'anosite', site3)) %>%
  mutate(site4 = gsub('nosite', 'anosite', site4)) %>%
  mutate(site5 = gsub('nosite', 'anosite', site5)) %>%
  mutate(site6 = gsub('nosite', 'anosite', site6))

ind_site_ind_back_fit <- ind_site_ind_back(subpool5_ncw)
ind_site_ind_back_sum <- tidy(ind_site_ind_back_fit) %>%
  filter(str_detect(term, '^site')) %>%
  mutate(term = gsub('consensus', '_consensus', term)) %>%
  mutate(term = gsub('weak', '_weak', term)) %>%
  separate(term, into = c('site', 'type'), sep = "_")

ind_site_ind_back_anova <- tidy(anova(ind_site_ind_back_fit)) %>%
  mutate(term_fctr = factor(term, levels = term))
ind_site_ind_back_p_r <- pred_resid(subpool5_ncw, ind_site_ind_back_fit)

p_ind_site_ind_back <- ggplot(ind_site_ind_back_p_r, 
                              aes(ave_med_ratio_norm, pred, fill = consensus)) +
  geom_point(shape = 21, alpha = 0.5) +
  scale_fill_viridis() +
  xlab('Measured expression') + ylab('Predicted expression') +
  annotate("text", x = 200, y = 0, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r$pred,
                                   ind_site_ind_back_p_r$ave_med_ratio_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_ind_site_ind_back_sum <- ggplot(ind_site_ind_back_sum, 
                                  aes(site, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_x_discrete(position = 'bottom') + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.ticks.x = element_blank()) +
  ylab('Weight')

p_ind_site_ind_back_anova <- ind_site_ind_back_anova %>%
  ggplot(aes(term_fctr, sumsq)) + 
  geom_bar(stat = 'identity') + 
  ylab('Sum of squares') +
  xlab('Model term') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())

p_ind_sit_ind_back_grid <- plot_grid(p_ind_site_ind_back, 
                                     p_ind_site_ind_back_sum, 
                                     p_ind_site_ind_back_anova,
                                     nrow = 3)

save_plot('plots/p_ind_sit_ind_back_grid.png', p_ind_sit_ind_back_grid,
          scale = 1.3, base_width = 5, base_height = 8)


#Compare sum and med analyses---------------------------------------------------

med_sum_join <- function(sum, med) {
  sum_select <- sum %>%
    select(name, ratio_br1, ratio_br2)
  med_select <- med %>%
    select(name, med_ratio_br1, med_ratio_br2)
  sum_med <- inner_join(sum_select, med_select, by = 'name') %>%
    select(-name)
  return(sum_med)
}

med_vs_sum <- med_sum_join(rep_1_2_sum, rep_1_2) %>%
  var_log10()

p_ave_med_vs_sum <- med_vs_sum %>%
  mutate(ave_sum = (ratio_br1 + ratio_br2)/2) %>%
  mutate(ave_med = (med_ratio_br1 + med_ratio_br2)/2) %>%
  ggplot(aes(ave_sum, ave_med)) +
  geom_point(alpha = 0.2)

my_points <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.2) +
    scale_x_continuous(limits = c(-3, 1.5), breaks = c(-3:1)) + 
    scale_y_continuous(limits = c(-3, 1.5), breaks = c(-3:1)) +
    annotation_logticks(sides = 'bl')
}

my_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(kernel = 'gaussian') +
    scale_x_continuous(limits = c(-3, 1.5), breaks = c(-3:1)) +
    scale_y_continuous(limits = c(-0.1, 1)) +
    annotation_logticks(sides = 'b')
}

p_med_vs_sum <- ggpairs(med_vs_sum, 
                        columnLabels = c('Sum Exp.\nRep. 1', 'Sum Exp.\nRep. 2',
                                         'Med Exp.\nRep. 1', 'Med Exp.\nRep. 2'),
                        lower = list(continuous = my_points),
                        diag = list(continuous = my_density)) +
  panel_border() + 
  theme(panel.grid.major = element_blank())

save_plot('plots/p_med_vs_sum.png', p_med_vs_sum, scale = 1.5)


#Do higher DNA BC reads contribute to high sum/median ratios?

#combine sum and median values, take sum/med

sum_over_med_join <- function(sum, med) {
  sum_select_1 <- sum %>%
    select(name, ratio_br1)
  med_select_1 <- med %>%
    select(name, med_ratio_br1)
  sum_med_1 <- inner_join(sum_select_1, med_select_1, by = 'name') %>%
    mutate(sum_over_med = ratio_br1/med_ratio_br1) %>%
    mutate(replicate = 1) %>%
    select(name, sum_over_med, replicate)
  sum_select_2 <- sum %>%
    select(name, ratio_br2)
  med_select_2 <- med %>%
    select(name, med_ratio_br2)
  sum_med_2 <- inner_join(sum_select_2, med_select_2, by = 'name') %>%
    mutate(sum_over_med = ratio_br2/med_ratio_br2) %>%
    mutate(replicate = 2) %>%
    select(name, sum_over_med, replicate)
  sum_med_rep <- rbind(sum_med_1, sum_med_2)
  return(sum_med_rep)
}

sum_over_med <- sum_over_med_join(rep_1_2_sum, rep_1_2)

#Combine sum/med with average DNA BC reads

DNA_norm_join_som <- function(df1, df2, df3) {
  DNA1 <- df1 %>%
    select(barcode, name, ave_norm_DNA) %>%
    mutate(replicate = 1)
  DNA2 <- df2 %>%
    select(barcode, name, ave_norm_DNA) %>%
    mutate(replicate = 2)
  DNA_rep <- rbind(DNA1, DNA2)
  DNA_join_som <- left_join(df3, DNA_rep, by = c('name', 'replicate'))
  return(DNA_join_som)
}

DNA_sum_over_med <- DNA_norm_join_som(bc_ave_DNA_RNA_1, bc_ave_DNA_RNA_2,
                                      sum_over_med)

p_DNA_sum_over_med <- ggplot(DNA_sum_over_med, aes(sum_over_med, ave_norm_DNA)) +
  facet_grid(~ replicate) +
  geom_point(alpha = 0.2) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  panel_border() + 
  xlab('Variant sum/med') +
  ylab('Average normalized DNA BC reads')

save_plot('plots/DNA_sum_over_med.png', p_DNA_sum_over_med, scale = 1.3)


#Combine sum/med with BC # per variant

BC_join_som <- function(df1, df2, df3) {
  BC1 <- df1 %>%
    select(name, barcodes_DNA) %>%
    mutate(replicate = 1)
  BC2 <- df2 %>%
    select(name, barcodes_DNA) %>%
    mutate(replicate = 2)
  BC_rep <- rbind(BC1, BC2)
  BC_join_som <- left_join(df3, BC_rep, by = c('name', 'replicate'))
  return(BC_join_som)
}

BC_sum_over_med <- BC_join_som(med_ratio_1, med_ratio_2, sum_over_med)

p_BC_sum_over_med <- ggplot(BC_sum_over_med, aes(sum_over_med, barcodes_DNA)) +
  facet_grid(~ replicate) +
  geom_point(alpha = 0.2) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  panel_border() + 
  xlab('Variant sum/med') +
  ylab('DNA Barcodes')

save_plot('plots/BC_sum_over_med.png', p_BC_sum_over_med, scale = 1.3)


#Finding variants for extra analyses--------------------------------------------

#Finding variants for Rishi

p_var_med_ratio <- ggplot(NULL, aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(data = filter(log2_rep_1_2), alpha = 0.2) +
  geom_point(data = filter(log2_rep_1_2, 
                           name == 'subpool5_weak_consensus_consensus_consensus_consensus_consensus_Vista Chr5:88673410-88674494'),
             fill = '#404788FF', color = 'black', 
             stroke = 0.2, size = 3, shape = 21) +
  geom_point(data = filter(log2_rep_1_2, 
                           name == 'subpool5_consensus_weak_consensus_weak_consensus_consensus_Vista Chr5:88673410-88674494' | name == 'subpool5_consensus_weak_consensus_no_site_consensus_consensus_Vista Chr5:88673410-88674494'),
             fill = '#1F968BFF', color = 'black', 
             stroke = 0.2, size = 3, shape = 21) +
  geom_point(data = filter(log2_rep_1_2, 
                           name == 'pGL4.29 Promega 1-63 + 1-87'),
             fill = '#95D840FF', color = 'black', 
             stroke = 0.2, size = 3, shape = 21) +
  annotation_logticks(scaled = TRUE) +
  xlab("log2 variant median\nnorm. RNA/DNA BR 1") +
  ylab("log2 variant median\nnorm. RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4), 
                     limits = c(-7.5, 4)) + 
  scale_y_continuous(breaks = c(-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4), 
                     limits = c(-7.5, 4)) + 
  annotate("text", x = -5, y = 2,
           label = paste(
             'r =', round(
               cor(
                 rep_1_2$med_ratio_br1, rep_1_2$med_ratio_br2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

high_br1 <- arrange(rep_1_2, desc(med_ratio_br1)) %>%
  head(1) %>%
  ungroup %>%
  select(name, most_common, med_ratio_br1, med_ratio_br2)

high_br2 <- arrange(rep_1_2, desc(med_ratio_br2)) %>%
  head(1) %>%
  ungroup %>%
  select(name, most_common, med_ratio_br1, med_ratio_br2)

high_br <- rbind(high_br1, high_br2) %>%
  write.table(
    "high_int_br.txt", 
    sep = '\t', row.names = FALSE)


#Random code snippets------------------------------------------------------------------------

sum(is.nan(subpool3_log10_norm$med_ratio_br1))

mutate(cv_med_ratio_norm = (
  sqrt(((med_ratio_br1_norm - ave_med_ratio_norm)^2)/(2-1)))/ave_med_ratio_norm) %>%
  filter(cv_med_ratio_norm <= 0.5)



