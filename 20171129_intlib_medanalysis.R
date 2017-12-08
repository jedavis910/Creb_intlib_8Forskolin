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


#Join DNA BC reads > 0 take average norm. reads and join RNA BC reads > 0 per BR--------------

ave_dna_join_rna_rep <- function(df1, df2, df3) {
  filter_reads_1 <- filter(df1, num_reads > 0)
  filter_reads_2 <- filter(df2, num_reads > 0)
  DNA_join <- inner_join(filter_reads_1, filter_reads_2, 
                         by = c("barcode", "name", "subpool", "most_common"), 
                         suffix = c("_DNA_tr1", "_DNA_tr2")
                         ) %>%
    mutate(ave_DNA_norm = (normalized_DNA_tr1 + normalized_DNA_tr2)/2)
  filter_reads_RNA <- filter(df3, num_reads > 0)
  DNA_RNA_join <- inner_join(DNA_join, filter_reads_RNA,
                             by = c("barcode", "name", "subpool", "most_common")
                             ) %>%
    rename(num_reads_RNA = num_reads) %>%
    rename(RNA_norm = normalized)
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
             suffix = c('_1br', '_2br'))

p_var_med_ratio <- ggplot(NULL, aes(med_ratio_1br, med_ratio_2br)) +
  geom_point(data = rep_1_2, alpha = 0.1) +
  geom_point(data = filter(rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 variant median\nnorm. RNA/DNA BR 1") +
  ylab("log10 variant median\nnorm. RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  annotate("text", x = -2, y = 1,
           label = paste(
             'r =', round(
               cor(
                 rep_1_2$med_ratio_1br, rep_1_2$med_ratio_2br,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

save_plot('plots/p_var_med_ratio.png', p_var_med_ratio)


