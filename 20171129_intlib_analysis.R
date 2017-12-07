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


#Combine bc technical replicates for each DNA biological replicate---------------------------

bc_dna_biol_rep <- function(df1, df2) {
  filter_reads_1 <- filter(df1, num_reads > 0)
  filter_reads_2 <- filter(df2, num_reads > 0)
  biol_rep <- inner_join(filter_reads_1, filter_reads_2, 
                         by = c("barcode", "name", "subpool", "most_common"), 
                         suffix = c("_1", "_2")
                         )
  average_tr <- biol_rep %>%
    mutate(average_norm = (normalized_1 + normalized_2)/2) %>%
    mutate(stdev_norm = ((normalized_1 - average_norm)^2)/2)
  print('processed dfs in order of technical replicates (1, 2) as (x, y)')
  return(average_tr)
}

bc_join_biol_rep_DNA_1 <- bc_dna_biol_rep(bc_join_DNA_1_1, bc_join_DNA_1_2)
bc_join_biol_rep_DNA_2 <- bc_dna_biol_rep(bc_join_DNA_2_1, bc_join_DNA_2_2)


#Determine variant counts by summing--------------------------------------------------------

#sum unique barcodes and normalized bc reads per variant. Output is total barcodes and sum 
#normalized reads per variant.

var_sum_bc_num_rna <- function(df1) {
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

var_sum_bc_num_dna <- function(df1) {
  bc_count <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = average_norm) %>%
    rename(sum = n)
  bc_sum <- inner_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")
  ) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_DNA_1 <- var_sum_bc_num_dna(bc_join_biol_rep_DNA_1)
variant_counts_DNA_2 <- var_sum_bc_num_dna(bc_join_biol_rep_DNA_2)
variant_counts_RNA_1 <- var_sum_bc_num_rna(bc_join_RNA_1)
variant_counts_RNA_2 <- var_sum_bc_num_rna(bc_join_RNA_2)


#Join RNA to DNA and determine expression from summing--------------------------------------

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets and determining 
#RNA/DNA per variant. Ratio is summed normalized reads of RNA over DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_RNA", "_DNA")
  ) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as RNA, y defined as DNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_1 <- var_expression(variant_counts_RNA_1, variant_counts_DNA_1)
RNA_DNA_2 <- var_expression(variant_counts_RNA_2, variant_counts_DNA_2)


#Combine BR and take log ratio of expression-------------------------------------------------

var_log_expression_rep <- function(df1, df2) {
  br <- inner_join(df1, df2,
                   by = c("name", "subpool", "most_common"), 
                   suffix = c("_1", "_2")
  )
  log_ratio_df <- br %>% 
    mutate_if(is.double, 
              funs(log10(.))
    )
  print('processed in format (br1, br2) in (x,y)')
  return(log_ratio_df)
}

log_rep_1_2 <- var_log_expression_rep(RNA_DNA_1, RNA_DNA_2)


#Replicate plot of summed expression between BR

p_br_var_expression <- ggplot(log_rep_1_2, aes(ratio_1, ratio_2)) +
  geom_point(data = log_rep_1_2, alpha = 0.3) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expr. Biol. Rep. 1") +
  ylab("Variant log10 Expr. Biol. Rep. 2") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  scale_y_continuous(breaks = c(-2, -1, 0, 1), limits = c(-2.5, 1.5)) + 
  annotate("text", x = -2, y = 1,
           label = paste(
             'r =', round(
               cor(
                 log_rep_1_2$ratio_1,log_rep_1_2$ratio_2,
                 use = "pairwise.complete.obs", method = "pearson"
               ), 2
             )
           )
  )

save_plot('plots/p_br_var_expression.png', p_br_var_expression, 
          scale = 0.9, base_height = 4, base_width = 4)


