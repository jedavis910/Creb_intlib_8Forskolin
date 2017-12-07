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


#Combine bc technical replicates for each DNA biological replicate and plot------------------

bc_dna_biol_rep <- function(df1, df2) {
  filter_reads_1 <- filter(df1, num_reads > 0)
  filter_reads_2 <- filter(df2, num_reads > 0)
  biol_rep <- inner_join(filter_reads_1, filter_reads_2, 
                       by = c("barcode", "name", "subpool", "most_common"), 
                       suffix = c("_1", "_2")
                       ) %>% 
    mutate_if(is.double, 
              funs(log10(.))
              )
  print('processed dfs in order of technical replicates (1, 2) as (x, y)')
  return(biol_rep)
}

bc_join_biol_rep_DNA_1 <- bc_dna_biol_rep(bc_join_DNA_1_1, bc_join_DNA_1_2)
bc_join_biol_rep_DNA_2 <- bc_dna_biol_rep(bc_join_DNA_2_1, bc_join_DNA_2_2)


#plot replicates for BC's similarly

p_bc_rep_DNA1 <- ggplot(bc_join_biol_rep_DNA_1, aes(normalized_1, normalized_2)) +
  geom_point(alpha = 0.1) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 Norm. BC reads Rep. 1") +
  ylab("Log10 Norm. BC reads Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3), limits = c(-1, 3)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3), limits = c(-1, 3)) + 
  annotate("text", x = 0, y = 2,
           label = paste(
             'r =', round(
               cor(
                 bc_join_biol_rep_DNA_1$normalized_1,
                 bc_join_biol_rep_DNA_1$normalized_2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_rep_DNA2 <- ggplot(bc_join_biol_rep_DNA_2, aes(normalized_1, normalized_2)) +
  geom_point(alpha = 0.1) +
  annotation_logticks(scaled = TRUE) +
  xlab("Log10 Norm. BC reads Rep. 1") +
  ylab("Log10 Norm. BC reads Rep. 2") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3), limits = c(-1, 3)) + 
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3), limits = c(-1, 3)) + 
  annotate("text", x = 0, y = 2,
           label = paste(
             'r =', round(
               cor(
                 bc_join_biol_rep_DNA_2$normalized_1,
                 bc_join_biol_rep_DNA_2$normalized_2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

p_bc_rep_grid <- plot_grid(p_bc_rep_DNA1, p_bc_rep_DNA2,
  labels = c("Biological Replicate 1", "Biological Replicate 2"),
  nrow = 1, ncol = 2, align = 'hv', hjust = -0.75, vjust = -0.25, scale = 0.9)

save_plot('plots/p_bc_rep_grid.png', 
          p_bc_rep_grid, base_height = 4, base_width = 8)


#Average normalized expression between TR per BR and determine variant counts by summing----

bc_dna_average <- function(df1) {
  average_tr <- df1 %>%
    mutate(average_norm = (normalized_1 + normalized_2)/2) %>%
    mutate(stdev_norm = ((normalized_1 - average_norm)^2)/2)
}

bc_ave_biol_rep_DNA_1 <- bc_dna_average(bc_join_biol_rep_DNA_1)
bc_ave_biol_rep_DNA_2 <- bc_dna_average(bc_join_biol_rep_DNA_2)

#sum unique barcodes and normalized bc reads per variant. Output is total barcodes and sum 
#normalized reads per variant. There is only 2 controls that are not represented by a BC read 
#in some of the samples, thus it is not present in joined tables later

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












