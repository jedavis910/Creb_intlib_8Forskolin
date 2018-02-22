library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)

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
#mapping andreplacing na with 0 reads.

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


#Join DNA and RNA BC------------------------------------------------------------

#Join DNA BC reads > 2, take average norm. reads and join to RNA. Take ratio of 
#RNA/DNA norm reads

ave_dna_join_rna_rep <- function(df1, df2, df3) {
  filter_reads_1 <- filter(df1, num_reads > 2)
  filter_reads_2 <- filter(df2, num_reads > 2)
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


#Barcodes per variant and median RNA/DNA----------------------------------------

#Count barcodes per variant per DNA and RNA, set minimum of 7 BC's per variant 
#in DNA sample, take median RNA/DNA per variant, find absolute deviation for
#each BC per variant then per variant determine the median absolute deviation.

ratio_bc_med_var <- function(df1) {
  bc_count_DNA <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df1 %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- inner_join(bc_count_DNA, bc_count_RNA, 
                           by = c('subpool', 'name', 'most_common'))
  med_ratio <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio))
  mad_ratio <- inner_join(df1, med_ratio, 
                          by = c('subpool', 'name', 'most_common')) %>%
    mutate(absdev = abs(ratio - med_ratio)) %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = median(absdev))
  med_mad <- inner_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio)) %>%
    mutate(mad_over_med = if_else(
      is.na(mad_over_med),
      as.double(0), 
      mad_over_med))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  return(bc_med)
}

med_ratio_1 <- ratio_bc_med_var(bc_ave_DNA_RNA_1)
med_ratio_2 <- ratio_bc_med_var(bc_ave_DNA_RNA_2)


#combine biological replicates, normalize to background---------------------------------------

#After combining, rename backgrounds to simplified names, make background column 
#(excluding controls), separate out background values in each dataset and left 
#join to original dataset. Normalize expression of each variant to its 
#background in that biological replicate. Determine average normalized 
#expression across biological replicates.

rep_1_2 <- inner_join(med_ratio_1, med_ratio_2,
             by = c("name", "subpool", "most_common"),
             suffix = c('_br1', '_br2'))

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


#determine the log(RNA/DNA) for each sample log2 is useful for replicate plots for expression 
#and log10 is useful for barcode read analysis. This is useful for replicate plots, but 
#further manipulations should use rep_1_2

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

int_back_norm_rep_1_2_log10 <- var_log10(int_back_norm_rep_1_2)

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
  filter(int_back_norm_rep_1_2_log10, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", 
    "fluff"), sep = "_"
  ) %>%
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

output_int <- int_back_norm_rep_1_2 %>%
  write.table(
    "int_back_norm_rep_1_2.txt", 
    sep = '\t', row.names = FALSE)

output_int_controls <- controls %>%
  write.table(
    'int_controls_rep_1_2.txt',
    sep = '\t', row.names = FALSE
  )


#plot variant statistics and rep plots-------------------------------------------------------

#comparing subpool expression

rep_1_2_0corr <- function(df1) {
  df1 <- df1 %>%
    mutate(med_ratio_br1 = if_else(
      med_ratio_br1 == as.double(0),
      as.double(0.01), 
      med_ratio_br1)) %>%
    mutate(med_ratio_br2 = if_else(
      med_ratio_br2 == as.double(0),
      as.double(0.01), 
      med_ratio_br2))
}

rep_1_2_0rm <- function(df1) {
  df1 <- df1 %>%
    filter(med_ratio_br1 != as.double(0)) %>%
    filter(med_ratio_br2 != as.double(0))
}

rep_1_2_log10 <- rep_1_2_0rm(rep_1_2) %>%
  var_log10()

test <- rep_1_2_0rm(rep_1_2) %>%
  filter(mad_over_med_br1 == 1 | mad_over_med_br2 == 1)

p_var_med_ratio <- ggplot(NULL, aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(data = filter(rep_1_2_log10, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(rep_1_2_log10, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = rep_1_2_log10, 
                 alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(rep_1_2_log10, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) + 
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant median\n RNA/DNA BR 1") +
  ylab("log10 variant median\n RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) + 
  scale_y_continuous(breaks = c(-3:1), limits = c(-3, 1.5))  +
  background_grid(major = 'xy', minor = 'none') + 
  annotate("text", x = -2, y = 0.5,
           label = paste('rtot =', 
                         round(cor(rep_1_2_log10$med_ratio_br1, 
                                   rep_1_2_log10$med_ratio_br2,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2))) +
  annotate("text", x = -2, y = 0, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(rep_1_2_log10, subpool == 'subpool5')$med_ratio_br1, 
               filter(rep_1_2_log10, subpool == 'subpool5')$med_ratio_br2,
               use = "pairwise.complete.obs", method = "pearson"), 2))) +
  annotate("text", x = -2, y = -0.5, color = '#2D708EFF',
         label = paste('r =', round(
           cor(
             filter(rep_1_2_log10, subpool == 'subpool3')$med_ratio_br1, 
             filter(rep_1_2_log10, subpool == 'subpool3')$med_ratio_br2,
             use = "pairwise.complete.obs", method = "pearson"), 2)))

save_plot('plots/p_var_med_ratio.png', p_var_med_ratio)

p_var_med_ratio_madmed1 <- ggplot(NULL, aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(data = filter(rep_1_2_log10, 
                           mad_over_med_br1 != as.double(0) & subpool != 'control'),
             color = 'black', alpha = 0.2) + 
  geom_point(data = filter(rep_1_2_log10, 
                           mad_over_med_br2 != as.double(0) & subpool != 'control'),
             color = 'black', alpha = 0.2) + 
  geom_point(data = filter(rep_1_2_log10, 
                         mad_over_med_br1 == as.double(0) & subpool != 'control'),
           color = 'red', alpha = 0.2) + 
  geom_point(data = filter(rep_1_2_log10, 
                           mad_over_med_br2 == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.2) + 
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


#Plot subpool expression features-----------------------------------------------------------

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

#BC num reads

p_BC_num_reads_viol_full <- ggplot(NULL) +
  geom_violin(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR1', y = num_reads_DNA_tr1, color = subpool)) +
  geom_violin(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR2', y = num_reads_DNA_tr2, color = subpool)) +
  geom_violin(data = bc_ave_DNA_RNA_1, 
              aes(x = 'RNA BR1', y = num_reads_RNA, color = subpool)) +
  geom_violin(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR1', y = num_reads_DNA_tr1, color = subpool)) +
  geom_violin(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR2', y = num_reads_DNA_tr2, color = subpool)) +
  geom_violin(data = bc_ave_DNA_RNA_2, 
              aes(x = 'RNA BR2', y = num_reads_RNA, color = subpool)) +
  xlab("") +
  ylab("Reads per BC")
  
save_plot('plots/BC_num_reads_viol_full.png',
          p_BC_num_reads_viol_full, scale = 2.2)

p_BC_num_reads_box_zoom <- ggplot(NULL) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR1', y = num_reads_DNA_tr1, color = subpool)) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'DNA BR1 TR2', y = num_reads_DNA_tr2, color = subpool)) +
  geom_boxplot(data = bc_ave_DNA_RNA_1, 
              aes(x = 'RNA BR1', y = num_reads_RNA, color = subpool)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR1', y = num_reads_DNA_tr1, color = subpool)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'DNA BR2 TR2', y = num_reads_DNA_tr2, color = subpool)) +
  geom_boxplot(data = bc_ave_DNA_RNA_2, 
              aes(x = 'RNA BR2', y = num_reads_RNA, color = subpool)) +
  xlab("") +
  ylab("Reads per BC") + ylim(0, 25)

save_plot('plots/BC_num_reads_box_zoom.png',
          p_BC_num_reads_box_zoom, scale = 2.2)

#BC normalized reads rep (need to log transform normalized reads first to compare, RNA and 
#ave DNA comparison required merging BC's that are present in each br (DNA and RNA))

norm_DNA_RNA_0corr <- function(df1) {
  df1 <- df1 %>%
    mutate(norm_DNA_tr1 = if_else(
      norm_DNA_tr1 == as.double(0),
      as.double(0.01), 
      norm_DNA_tr1)) %>%
    mutate(norm_DNA_tr2 = if_else(
      norm_DNA_tr2 == as.double(0),
      as.double(0.01), 
      norm_DNA_tr2)) %>%
    mutate(norm_RNA = if_else(
      norm_RNA == as.double(0),
      as.double(0.01), 
      norm_RNA))
}

bc_log_ave_DNA_RNA_1 <- norm_DNA_RNA_0corr(bc_ave_DNA_RNA_1) %>%
  var_log10()
bc_log_ave_DNA_RNA_2 <- norm_DNA_RNA_0corr(bc_ave_DNA_RNA_2) %>%
  var_log10()
bc_log_ave_br <- inner_join(bc_log_ave_DNA_RNA_1, bc_log_ave_DNA_RNA_2, 
                            by = c("barcode", "name", "subpool", "most_common"),
                            suffix = c('_br1', '_br2'))

p_bc_rep_DNA1 <- ggplot(data = NULL, aes(norm_DNA_tr1, norm_DNA_tr2)) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_1, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_1, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_1, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads TR 1") +
  ylab("Log10 norm. BC reads TR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  scale_y_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  annotate("text", x = -2, y = 3,
           label = paste(
             'r =', round(
               cor(bc_log_ave_DNA_RNA_1$norm_DNA_tr1, 
                   bc_log_ave_DNA_RNA_1$norm_DNA_tr2,
                   use = "pairwise.complete.obs", method = "pearson"), 2)))

p_bc_rep_DNA2 <- ggplot(data = NULL, aes(norm_DNA_tr1, norm_DNA_tr2)) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_2, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_2, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_DNA_RNA_2, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads TR 1") +
  ylab("Log10 norm. BC reads TR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  scale_y_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  annotate("text", x = -2, y = 3,
           label = paste(
             'r =', round(
               cor(bc_log_ave_DNA_RNA_2$norm_DNA_tr1, 
                   bc_log_ave_DNA_RNA_2$norm_DNA_tr2,
                   use = "pairwise.complete.obs", method = "pearson"), 2)))

p_bc_rep_ave <- ggplot(data = NULL, aes(ave_norm_DNA_br1, ave_norm_DNA_br2)) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 norm. BC reads BR 1") +
  ylab("Ave log10 norm. BC reads BR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  scale_y_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  annotate("text", x = -2, y = 3,
           label = paste('r =', round(cor(bc_log_ave_br$ave_norm_DNA_br1,
                                          bc_log_ave_br$ave_norm_DNA_br2,
                                          use = "pairwise.complete.obs", 
                                          method = "pearson"), 2)))

p_bc_rep_RNA <- ggplot(data = NULL, aes(norm_RNA_br1, norm_RNA_br2)) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads BR 1") +
  ylab("Log10 norm. BC reads BR 2") +
  background_grid(major = 'xy', minor = 'none') + 
  scale_x_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  scale_y_continuous(breaks = c(-2:3), limits = c(-2.5, 3.5)) + 
  annotate("text", x = -2, y = 3, 
           label = paste('r =', round(cor(bc_log_ave_br$norm_RNA_br1,
                                          bc_log_ave_br$norm_RNA_br2,
                                          use = "pairwise.complete.obs", 
                                          method = "pearson"), 2)))

p_bc_rep_grid <- plot_grid(p_bc_rep_DNA1, p_bc_rep_DNA2, p_bc_rep_ave, 
                           p_bc_rep_RNA,
                           labels = c(
                             'BR1 DNA', 'BR2 DNA', 'Ave DNA', '      RNA'),
                           nrow = 2, ncol = 2, align = 'hv', 
                           hjust = -2, vjust = 0.5, scale = 0.9)

save_plot('plots/p_bc_rep_grid.png', p_bc_rep_grid, 
          base_height = 7, base_width = 10)


#Sum analysis-------------------------------------------------------------------

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = norm) %>%
    rename(sum = n)
  bc_sum <- inner_join(variant_sum, bc_count, 
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


ave_dna_join_rna_var <- function(df1, df2, df3) {
  DNA_join <- inner_join(df1, df2, 
                         by = c("name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")) %>%
    mutate(ave_sum_DNA = (sum_DNA_tr1 + sum_DNA_tr2)/2) %>%
    filter(barcodes_DNA_tr1 > 7) %>%
    filter(barcodes_DNA_tr2 > 7)
  DNA_RNA_join <- left_join(DNA_join, df3,
                            by = c("name", "subpool", "most_common")) %>%
    rename(sum_RNA = sum) %>%
    rename(barcodes_RNA = barcodes) %>%
    mutate(ratio = sum_RNA/ave_sum_DNA) %>%
    mutate(sum_RNA = if_else(is.na(sum_RNA),
                             0,
                             sum_RNA)) %>%
    mutate(ratio = if_else(is.na(ratio), 
                          0, 
                          ratio)) %>%
    mutate(barcodes_RNA = if_else(is.na(barcodes_RNA),
                                  as.integer(0),
                                  barcodes_RNA))
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

rep_1_2_sum_0corr <- function(df1) {
  df1 <- df1 %>%
    mutate(ratio_br1 = if_else(
      ratio_br1 == as.double(0),
      as.double(0.01), 
      ratio_br1)) %>%
    mutate(ratio_br2 = if_else(
      ratio_br2 == as.double(0),
      as.double(0.01), 
      ratio_br2))
}

rep_1_2_sum_0rm <- function(df1) {
  df1 <- df1 %>%
    filter(ratio_br1 != as.double(0)) %>%
    filter(ratio_br2 != as.double(0))
}

rep_1_2_sum_log10 <- rep_1_2_sum_0rm(rep_1_2_sum) %>%
  var_log10()

p_var_sum_ratio <- ggplot(NULL, aes(ratio_br1, ratio_br2)) +
  geom_point(data = filter(rep_1_2_sum_log10, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(rep_1_2_sum_log10, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = rep_1_2_sum_log10, 
                 alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(rep_1_2_sum_log10, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) + 
  annotation_logticks(scaled = TRUE) +
  background_grid(major = 'xy', minor = 'none') + 
  xlab("log10 variant sum\n RNA/DNA BR 1") +
  ylab("log10 variant sum\n RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) + 
  scale_y_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) +
  annotate("text", x = -2, y = 0.5,
           label = paste('rtot =', 
                         round(cor(rep_1_2_sum_log10$ratio_br1, 
                                   rep_1_2_sum_log10$ratio_br2,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2))) +
  annotate("text", x = -2, y = 0, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(rep_1_2_sum_log10, subpool == 'subpool5')$ratio_br1, 
               filter(rep_1_2_sum_log10, subpool == 'subpool5')$ratio_br2,
               use = "pairwise.complete.obs", method = "pearson"), 2))) +
  annotate("text", x = -2, y = -0.5, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(rep_1_2_sum_log10, subpool == 'subpool3')$ratio_br1, 
               filter(rep_1_2_sum_log10, subpool == 'subpool3')$ratio_br2,
               use = "pairwise.complete.obs", method = "pearson"), 2)))

save_plot('plots/p_var_sum_ratio.png', p_var_sum_ratio)


#Compare sum and med analyses---------------------------------------------------

ave_rep_sum_0rm <- rep_1_2_sum_0rm(rep_1_2_sum) %>%
  mutate(ave_sum_ratio = (ratio_br1 + ratio_br2)/2) %>%
  var_log10()

ave_rep_med_0rm <- rep_1_2_0rm(rep_1_2) %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2) %>%
  var_log10()

rep_sum_med <- inner_join(ave_rep_sum_0rm, ave_rep_med_0rm, 
                          by = c('subpool', 'name', 'most_common'),
                          suffix = c('_sum', '_med'))

p_var_sum_med <- ggplot(NULL, aes(ave_sum_ratio, ave_med_ratio)) +
  geom_point(data = filter(rep_sum_med, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(rep_sum_med, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = rep_sum_med, 
                 alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(rep_sum_med, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) + 
  annotation_logticks(scaled = TRUE) +
  xlab("log10 averagevariant\nsum RNA/DNA ") +
  ylab("log10 average variant\nmedian RNA/DNA ") +
  scale_x_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) + 
  scale_y_continuous(breaks = c(-3:1), limits = c(-3, 1.5))  +
  background_grid(major = 'xy', minor = 'none') + 
  annotate("text", x = -2, y = 0.5,
           label = paste('rtot =', 
                         round(cor(rep_sum_med$ave_sum_ratio, 
                                   rep_sum_med$ave_med_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2))) +
  annotate("text", x = -2, y = 0, color = '#482677FF',
           label = paste('r =', round(
             cor(
               filter(rep_sum_med, subpool == 'subpool5')$ave_sum_ratio, 
               filter(rep_sum_med, subpool == 'subpool5')$ave_med_ratio,
               use = "pairwise.complete.obs", method = "pearson"), 2))) +
  annotate("text", x = -2, y = -0.5, color = '#2D708EFF',
           label = paste('r =', round(
             cor(
               filter(rep_sum_med, subpool == 'subpool3')$ave_sum_ratio, 
               filter(rep_sum_med, subpool == 'subpool3')$ave_med_ratio,
               use = "pairwise.complete.obs", method = "pearson"), 2)))

save_plot('plots/p_var_sum_med.png', p_var_sum_med)

p_var_sum_med_madmed1 <- ggplot(NULL, aes(ave_sum_ratio, ave_med_ratio)) +
  geom_point(data = filter(rep_sum_med, 
                           mad_over_med_br1 != as.double(0) & subpool != 'control'),
             color = 'black', alpha = 0.2) + 
  geom_point(data = filter(rep_sum_med, 
                           mad_over_med_br2 != as.double(0) & subpool != 'control'),
             color = 'black', alpha = 0.2) + 
  geom_point(data = filter(rep_sum_med, 
                           mad_over_med_br1 == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.2) + 
  geom_point(data = filter(rep_sum_med, 
                           mad_over_med_br2 == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.2) + 
  annotation_logticks(scaled = TRUE) +
  xlab("log10 averagevariant\nsum RNA/DNA ") +
  ylab("log10 average variant\nmedian RNA/DNA ") +
  scale_x_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) + 
  scale_y_continuous(breaks = c(-3:1), limits = c(-3, 1.5)) +
  annotate("text", x = -2, y = 0.5,
           label = paste('rtot =', 
                         round(cor(rep_sum_med$ave_sum_ratio, 
                                   rep_sum_med$ave_med_ratio,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))


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



