library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)

#Written for the analysis of integrated library expression at 8 µM. 2 biological replicates
#and 2 technical replicates each for DNA amplification
#D##_BC: DNA biological replicate # technical replicate #
#R#8_BC: RNA biological replicate # at 8 µM


#Load index and bcmap files------------------------------------------------------------------

bc_DNA_1_1 <- read_tsv('BCreads_txts/D11_BC.txt')
bc_DNA_1_2 <- read_tsv('BCreads_txts/D12_BC.txt')
bc_DNA_2_1 <- read_tsv('BCreads_txts/D21_BC.txt')
bc_DNA_2_2 <- read_tsv('BCreads_txts/D22_BC.txt')
bc_RNA_1 <- read_tsv('BCreads_txts/R18_BC.txt')
bc_RNA_2 <- read_tsv('BCreads_txts/R28_BC.txt')


#Load barcode mapping table, remember sequences are rcomp due to sequencing format

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'
                        ), 
                        skip = 1) %>%
  select(-fluff)


#pick out SP3 and SP5 in the bcmap that were used in this assay

SP3_SP5_map <- barcode_map %>%
  mutate(subpool = ifelse(
    startsWith(name, 'subpool'), 
    substr(name, 1, 8), 
    'control'
  )
  ) %>%
  filter(subpool != 'subpool2') %>%
  filter(subpool != 'subpool4')


#Join reads to bcmap------------------------------------------------------------------------

#Join BC reads to BC mapping, keeping the reads only appearing in barcode mapping and 
#replacing na with 0 reads.

bc_map_join_bc <- function(df1, df2) {
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(num_reads = if_else(
      is.na(num_reads), 
      as.integer(0), 
      num_reads
    )
    ) %>%
    mutate(norm = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  return(keep_bc)
}

bc_join_DNA_1_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_1)
bc_join_DNA_1_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_2)
bc_join_DNA_2_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_1)
bc_join_DNA_2_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_2)
bc_join_RNA_1 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_1)
bc_join_RNA_2 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_2)


#Join DNA BC reads > 0 take average norm. reads and join RNA BC reads > 0 per BR--------------

ave_dna_join_rna_rep <- function(df1, df2, df3) {
  filter_reads_1 <- filter(df1, num_reads > 0)
  filter_reads_2 <- filter(df2, num_reads > 0)
  DNA_join <- inner_join(filter_reads_1, filter_reads_2, 
                         by = c("barcode", "name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")
                         ) %>%
    mutate(ave_DNA_norm = (norm_DNA_tr1 + norm_DNA_tr2)/2)
  filter_reads_RNA <- filter(df3, num_reads > 0)
  DNA_RNA_join <- inner_join(DNA_join, filter_reads_RNA,
                             by = c("barcode", "name", "subpool", "most_common")
                             ) %>%
    rename(num_reads_RNA = num_reads) %>%
    rename(RNA_norm = norm)
  print('processed dfs in order of (DNA tr1, DNA tr2, RNA) in bc_dna_biol_rep(df1, df2, df3)')
  return(DNA_RNA_join)
}

bc_ave_DNA_RNA_1 <- ave_dna_join_rna_rep(bc_join_DNA_1_1, bc_join_DNA_1_2, bc_join_RNA_1)
bc_ave_DNA_RNA_2 <- ave_dna_join_rna_rep(bc_join_DNA_2_1, bc_join_DNA_2_2, bc_join_RNA_2)


#Determine RNA/DNA per BC in each set, determine median RNA/DNA per variant------------------

ratio_bc_med_var <- function(df1) {
  bc_count <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  med_ratio <- df1 %>%
    mutate(ratio = RNA_norm/ave_DNA_norm) %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio)) %>%
    mutate(med_ratio = log10(med_ratio))
  bc_med <- inner_join(med_ratio, bc_count, 
                       by = c("name", "subpool", "most_common")
  )
  return(bc_med)
}

med_ratio_1 <- ratio_bc_med_var(bc_ave_DNA_RNA_1)
med_ratio_2 <- ratio_bc_med_var(bc_ave_DNA_RNA_2)


#combine biological replicates---------------------------------------------------------------

rep_1_2 <- inner_join(med_ratio_1, med_ratio_2,
             by = c("name", "subpool", "most_common"),
             suffix = c('_br1', '_br2'))


#Separate into subpools----------------------------------------------------------------------

#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) that vary in 
#distance from one another by 0 (no inner flanks), 5, 10, 15, 20 and 70 bp (all but 0 appear 
#as -4 bp spacing). Each site distance combination is then moved along the backgrounds at 1 
#bp increments starting from closest to the minP. Separation lists the spacing between sites,
#distance (start of consensus and flanks) and the background. Added 2 to all distances to 
#measure to start of BS and not to flank. Added 4 to all spacing but 0 to measure difference 
#between start of sites. Also took average of log10 med BC expression between biological 
#replicates for plotting

subpool3 <- 
  filter(rep_1_2, subpool == "subpool3") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, 
           into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "background"),
           sep = "_", convert = TRUE
           ) %>%
  select(-subpool, -fluff2, -fluff3) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)
         ) %>%
  mutate(dist = dist + 2) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing)) %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from furthest to 
#the minP. These sites are filled with sites of either the consensus site, a weak site or no 
#site. Both the weak and consensus sites are flanked by the same flanking sequence. Also took
#average of log10 med BC expression between biological replicates for plotting

subpool5 <- 
  filter(rep_1_2, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", 
    "background"), sep = "_"
  ) %>%
  select(-subpool) %>%
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
                  'weak', site_combo)) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(ave_med_ratio = (med_ratio_br1 + med_ratio_br2)/2)

controls <- 
  filter(rep_1_2, subpool == "control") %>%
  ungroup ()


#Plot subpool expression features-----------------------------------------------------------

#Subpool 3

p_subpool3_spa_back <- ggplot(subpool3) + 
  geom_point(aes(dist, med_ratio_br1), alpha = 0.7, size = 1.5, color = '#287D8EFF') +
  geom_point(aes(dist, med_ratio_br2), alpha = 0.7, size = 1.5, color = '#95D840FF') +
  facet_grid(spacing ~ background) + 
  geom_smooth(aes(dist, ave_med_ratio), 
              span = 0.1, size = 0.7, color = '#481567FF', se = FALSE
              ) +
  ylab('log10 median BC expression') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
    breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/p_subpool3_spa_back.png', p_subpool3_spa_back, 
          base_width = 46, base_height = 17, scale = 0.35)

#Subpool 5

p_subpool5_cons_mix_weak_sitenum <- ggplot(NULL, aes(as.factor(total_sites),ave_med_ratio)) + 
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

p_subpool5_cons_mix_sitenum <- ggplot(NULL, aes(site_type, ave_med_ratio)) +
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

bc_log_ave_DNA_RNA_1 <- bc_ave_DNA_RNA_1 %>%
  mutate_if(is.double, funs(log10(.)))

bc_log_ave_DNA_RNA_2 <- bc_ave_DNA_RNA_2 %>%
  mutate_if(is.double, funs(log10(.)))

bc_log_ave_br <- inner_join(bc_log_ave_DNA_RNA_1, bc_log_ave_DNA_RNA_2, 
                            by = c("barcode", "name", "subpool", "most_common"),
                            suffix = c('_br1', '_br2')
                            )

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
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = -1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 bc_log_ave_DNA_RNA_1$norm_DNA_tr1,
                 bc_log_ave_DNA_RNA_1$norm_DNA_tr2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

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
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = -1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 bc_log_ave_DNA_RNA_2$norm_DNA_tr1,
                 bc_log_ave_DNA_RNA_2$norm_DNA_tr2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )


p_bc_rep_ave <- ggplot(data = NULL, aes(ave_DNA_norm_br1, ave_DNA_norm_br2)) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Ave log10 norm. BC reads BR 1") +
  ylab("Ave log10 norm. BC reads BR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = -1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 bc_log_ave_br$ave_DNA_norm_br1,
                 bc_log_ave_br$ave_DNA_norm_br2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_rep_RNA <- ggplot(data = NULL, aes(RNA_norm_br1, RNA_norm_br2)) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_point(data = filter(bc_log_ave_br, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.2) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 norm. BC reads BR 1") +
  ylab("Log10 norm. BC reads BR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = -1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 bc_log_ave_br$RNA_norm_br1,
                 bc_log_ave_br$RNA_norm_br2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_rep_grid <- plot_grid(p_bc_rep_DNA1, p_bc_rep_DNA2, p_bc_rep_ave, p_bc_rep_RNA,
                           labels = c('BR1 DNA', 'BR2 DNA', 'Ave DNA', '      RNA'),
                           nrow = 2, ncol = 2, align = 'hv', 
                           hjust = -2, vjust = 0.5, scale = 0.9)

save_plot('plots/p_bc_rep_grid.png', p_bc_rep_grid, base_height = 7, base_width = 10)


#plot variant statistics and rep plots-------------------------------------------------------

#comparing subpool expression

p_var_med_ratio <- ggplot(NULL, aes(med_ratio_br1, med_ratio_br2)) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.2) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.2) +
  geom_density2d(data = rep_1_2, alpha = 0.7, color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(rep_1_2, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.5) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) + 
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant median\nnorm. RNA/DNA BR 1") +
  ylab("log10 variant median\nnorm. RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  annotate("text", x = -2, y = 1,
           label = paste(
             'r =', round(
               cor(
                 rep_1_2$med_ratio_br1, rep_1_2$med_ratio_br2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

save_plot('plots/p_var_med_ratio.png', p_var_med_ratio)

#comparing BC number to med expression

p_var_med_bc_1 <- ggplot(NULL, aes(barcodes_br1, med_ratio_br1)) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.3) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.3) +
  geom_point(data = filter(rep_1_2, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.6) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.8) + 
  annotation_logticks(scaled = TRUE, sides = 'l') +
  ylab("log10 variant median\nnorm. RNA/DNA") +
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  xlab("Barcodes per variant") +
  annotate("text", x = 25, y = -2,
           label = paste(
             'r =', round(
               cor(
                 rep_1_2$med_ratio_br1, rep_1_2$barcodes_br1,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_med_bc_2 <- ggplot(NULL, aes(barcodes_br2, med_ratio_br2)) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool5'),
             color = '#482677FF', alpha = 0.3) +
  geom_point(data = filter(rep_1_2, subpool == 'subpool3'),
             color = '#2D708EFF', alpha = 0.3) +
  geom_point(data = filter(rep_1_2, subpool == 'control'),
             color = '#3CBB75FF', alpha = 0.6) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.8) + 
  annotation_logticks(scaled = TRUE, sides = 'l') +
  ylab("log10 variant median\nnorm. RNA/DNA") +
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  xlab("Barcodes per variant") +
  annotate("text", x = 25, y = -2,
           label = paste(
             'r =', round(
               cor(
                 rep_1_2$med_ratio_br2, rep_1_2$barcodes_br2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_med_bc <- plot_grid(p_var_med_bc_1, p_var_med_bc_2,
                           labels = c('BR1', 'BR2'),
                           nrow = 1, ncol = 2, align = 'hv', 
                           hjust = -3, vjust = 0.5, scale = 0.9)

save_plot('plots/p_var_med_bc.png', p_var_med_bc, base_height = 4, base_width = 10)


