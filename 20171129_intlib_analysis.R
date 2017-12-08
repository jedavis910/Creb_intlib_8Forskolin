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
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  return(keep_bc)
}

bc_join_DNA_1_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_1)
bc_join_DNA_1_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_1_2)
bc_join_DNA_2_1 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_1)
bc_join_DNA_2_2 <- bc_map_join_bc(SP3_SP5_map, bc_DNA_2_2)
bc_join_RNA_1 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_1)
bc_join_RNA_2 <- bc_map_join_bc(SP3_SP5_map, bc_RNA_2)


#Determine variant counts by summing--------------------------------------------------------

#sum unique barcodes and normalized bc reads per variant. Output is total barcodes and sum 
#normalized reads per variant.

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = normalized) %>%
    rename(sum = n)
  bc_sum <- inner_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")
  ) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA_1_1 <- var_sum_bc_num(bc_join_DNA_1_1)
variant_counts_DNA_1_2 <- var_sum_bc_num(bc_join_DNA_1_2)
variant_counts_DNA_2_1 <- var_sum_bc_num(bc_join_DNA_2_1)
variant_counts_DNA_2_2 <- var_sum_bc_num(bc_join_DNA_2_2)
variant_counts_RNA_1 <- var_sum_bc_num(bc_join_RNA_1)
variant_counts_RNA_2 <- var_sum_bc_num(bc_join_RNA_2)


#Join DNA replicates together and take average of summed normalized reads-------------------

#combine DNA tr, only keeping instances in both and take average of summed normalized reads. 
#Combine DNA with RNA per br, only keeping instances in both, and determine expression as the
#ratio of summed normalized RNA reads over the average of summed normalized DNA reads

DNA_RNA_expression <- function(df1, df2, df3) {
  DNA_tr_join <- inner_join(df1, df2, 
                            by = c('subpool', 'name', 'most_common'),
                            suffix = c('_DNA_1tr', '_DNA_2tr')
                            ) %>%
    mutate(ave_DNA_sum = (sum_DNA_1tr + sum_DNA_2tr)/2)
  DNA_RNA_join_expression <- inner_join(DNA_tr_join, df3,
                                        by = c('subpool', 'name', 'most_common')
                                        ) %>%
    rename(sum_RNA = sum) %>%
    rename(barcodes_RNA = barcodes) %>%
    mutate(ratio = sum_RNA/ave_DNA_sum)
  print('processed each biological replicate in order of (DNA_tr_1, DNA_tr_2, RNA) in 
        DNA_RNA_expression(df1, df2, df3)')
  return(DNA_RNA_join_expression)
} 

RNA_DNA_1 <- DNA_RNA_expression(
  variant_counts_DNA_1_1, variant_counts_DNA_1_2, variant_counts_RNA_1
  )
RNA_DNA_2 <- DNA_RNA_expression(
  variant_counts_DNA_2_1, variant_counts_DNA_2_2, variant_counts_RNA_2
)


#Combine BR and take log ratio of expression-------------------------------------------------

var_log_expression_rep <- function(df1, df2) {
  br <- inner_join(df1, df2,
                   by = c("name", "subpool", "most_common"), 
                   suffix = c("_1br", "_2br")
  )
  log_ratio_df <- br %>% 
    mutate_if(is.double, 
              funs(log10(.))
    ) %>%
    mutate(sd_DNA_norm_1br = ((sum_DNA_1tr_1br - ave_DNA_sum_1br)^2)/2) %>%
    mutate(sd_DNA_norm_2br = ((sum_DNA_1tr_2br - ave_DNA_sum_2br)^2)/2) %>%
    mutate(ave_DNA_bc_1br = (barcodes_DNA_1tr_1br + barcodes_DNA_2tr_1br)/2) %>%
    mutate(ave_DNA_bc_2br = (barcodes_DNA_1tr_2br + barcodes_DNA_2tr_2br)/2)
  print('processed in format (br1, br2) in (x,y)')
  return(log_ratio_df)
}

log_rep_1_2 <- var_log_expression_rep(RNA_DNA_1, RNA_DNA_2)


#Replicate plots-----------------------------------------------------------------------------

#BC number correlation between tr and br

p_bc_DNA_1br <- ggplot(log_rep_1_2, aes(barcodes_DNA_1tr_1br, barcodes_DNA_2tr_1br)) +
  geom_point(alpha = 0.1) + 
  xlab("BC's per variant TR 1") +
  ylab("BC's per variant TR 2") +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  annotate("text", x = 25, y = 75,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$barcodes_DNA_1tr_1br,log_rep_1_2$barcodes_DNA_2tr_1br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_DNA_2br <- ggplot(log_rep_1_2, aes(barcodes_DNA_1tr_2br, barcodes_DNA_2tr_2br)) +
  geom_point(alpha = 0.1) + 
  xlab("BC's per variant TR 1") +
  ylab("BC's per variant TR 2") +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  annotate("text", x = 25, y = 75,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$barcodes_DNA_1tr_2br,log_rep_1_2$barcodes_DNA_2tr_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_DNA_ave <- ggplot(log_rep_1_2, aes(ave_DNA_bc_1br, ave_DNA_bc_2br)) +
  geom_point(alpha = 0.1) + 
  xlab("Ave BC's per variant BR 1") +
  ylab("Ave BC's per variant BR 2") +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  annotate("text", x = 25, y = 75,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$ave_DNA_bc_1br,log_rep_1_2$ave_DNA_bc_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_RNA <- ggplot(log_rep_1_2, aes(barcodes_RNA_1br, barcodes_RNA_2br)) +
  geom_point(alpha = 0.1) + 
  xlab("BC's per variant BR 1") +
  ylab("BC's per variant BR 2") +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100, 120, 140, 160), 
                     limits = c(0, 160)) + 
  annotate("text", x = 25, y = 75,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$barcodes_RNA_1br,log_rep_1_2$barcodes_RNA_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_tr_br <- plot_grid(p_bc_DNA_1br, p_bc_DNA_2br, p_bc_DNA_ave, p_bc_RNA,
                        labels = c('BR1 DNA', 'BR2 DNA', 'Ave DNA', '      RNA'),
                        nrow = 2, ncol = 2, align = 'hv', 
                        hjust = -2, vjust = 0.5, scale = 0.9)

save_plot('plots/p_bc_tr_br.png', p_bc_tr_br, base_height = 7, base_width = 8)


#Correlation of summed normalized reads

p_var_sum_DNA_1br <- ggplot(NULL, aes(sum_DNA_1tr_1br, sum_DNA_2tr_1br)) +
  geom_point(data = log_rep_1_2, alpha = 0.1) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant summed norm. reads\nTR 1") +
  ylab("log10 variant summed norm. reads\nTR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = 1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$sum_DNA_1tr_1br,log_rep_1_2$sum_DNA_2tr_1br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_sum_DNA_2br <- ggplot(NULL, aes(sum_DNA_1tr_2br, sum_DNA_2tr_2br)) +
  geom_point(data = log_rep_1_2, alpha = 0.1) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant summed norm. reads\nTR 1") +
  ylab("log10 variant summed norm. reads\nTR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = 1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$sum_DNA_1tr_2br,log_rep_1_2$sum_DNA_2tr_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_sum_DNA_ave <- ggplot(NULL, aes(ave_DNA_sum_1br, ave_DNA_sum_2br)) +
  geom_point(data = log_rep_1_2, alpha = 0.1) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 ave variant summed norm. reads\nBR 1") +
  ylab("log10 ave variant summed norm. reads\nBR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = 1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$ave_DNA_sum_1br,log_rep_1_2$ave_DNA_sum_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_sum_RNA <- ggplot(NULL, aes(sum_RNA_1br, sum_RNA_2br)) +
  geom_point(data = log_rep_1_2, alpha = 0.1) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant summed norm. reads\nBR 1") +
  ylab("log10 variant summed norm. reads\nBR 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1.5, 4)) + 
  annotate("text", x = 1, y = 3,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$sum_RNA_1br,log_rep_1_2$sum_RNA_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_sum_ratio <- ggplot(NULL, aes(ratio_1br, ratio_2br)) +
  geom_point(data = log_rep_1_2, alpha = 0.1) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant summed norm. reads\nRNA/DNA BR 1") +
  ylab("log10 variant summed norm. reads\nRNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  annotate("text", x = -1, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$ratio_1br,log_rep_1_2$ratio_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_var_sum_tr_br <- plot_grid(
  p_var_sum_DNA_1br, p_var_sum_DNA_2br, p_var_sum_DNA_ave, p_var_sum_RNA, p_var_sum_ratio,
  labels = c('BR1 DNA', 'BR2 DNA', ' Ave DNA', '     RNA', 'RNA/DNA'),
  nrow = 3, ncol = 2, align = 'hv', hjust = -1.5, vjust = 0.5, scale = 0.9
  )

save_plot('plots/p_var_sum_tr_br.png', p_var_sum_tr_br, base_height = 14, base_width = 10)



