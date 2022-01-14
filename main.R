subset_mat <- function(mat) {
  counts <- read.table(mat, header=TRUE, row.names='gene')
  subset_mat <- counts[, c('vP0_1', 'vP0_2', 'vAd_1', 'vAd_2')]
  return(subset_mat)
  
}

make_coldata_subset <- function(metadata, c1, c2) {
  coldata <- data.frame(readr::read_csv(metadata) %>% 
                          select(samplename, timepoint) %>% 
                          filter(timepoint %in% c(c1, c2)) %>% 
                          mutate_at(vars(timepoint), factor))
  
  return(coldata)
}

#results as a list to return multiple values in R
deseq2_run <- function(count_mat, meta_data, design_formula) {
  dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=meta_data, design = design_formula)
  dds$timepoint <- relevel(dds$timepoint, ref='vP0')
  dds <- DESeq(dds)

  res <- results(dds, contrast=c('timepoint', 'vAd', 'vP0'))
  res_ordered <- res[order(res$padj),]
  results <- c(dds_obj = dds, res_mat = res_ordered)
  
return(results)
}

label_results <- function(results) {
res_labeled <- results %>%
  as.data.frame() %>%
  as_tibble(rownames='genes') %>%
  mutate(volc_plot_status = case_when(log2FoldChange > 0 & padj < .10 ~ 'Up',
                                      log2FoldChange < 0 & padj < .10 ~ 'Down',
                                      TRUE ~ 'NS'))
return(res_labeled)
}

plot_pvals <- function(res_labeled) {
pval_plot <- res_labeled %>% 
  ggplot(aes(pvalue)) + 
  geom_histogram(bins=50, color='black', fill='lightblue') + 
  theme_linedraw() + 
  ggtitle('Histogram of raw pvalues obtained from DE analysis') + 
  theme(plot.title = element_text(hjust = 0.5))
return(pval_plot)
}

plot_log2fc <- function(res_labeled) {
logfc_plot <- res_labeled %>% 
  filter(padj < .1) %>% 
  ggplot(aes(log2FoldChange)) + 
  geom_histogram(bins=100, color='black', fill='light blue') + 
  theme_linedraw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle('Histogram of Log2FoldChanges for DE Genes')
return(logfc_plot)
}

plot_volcano <- function(res_labeled) {
volcano_plot <- res_labeled %>% 
  ggplot() + 
  geom_point(mapping=aes(x=log2FoldChange, y=-log10(padj), color=volc_plot_status)) + 
  geom_hline(yintercept = -log10(0.1), linetype = "dashed")  + 
  theme_linedraw()
return(volcano_plot)
}

return_top_ten <- function(res_labeled) {
top_ten <- res_labeled %>% 
  top_n(-10, padj) %>% 
  dplyr::select(genes)
return(top_ten)
}

deseq2_norm_counts <- function(deseq2_obj) {
norm_counts <- counts(deseq2_obj, normalized=T)
return(norm_counts)
}

scatter_norm_counts <- function(norm_counts, top_ten){
norm_counts_tb <- norm_counts %>% as_tibble(rownames='genes')
scatter_norm <- top_ten %>% 
  left_join(norm_counts_tb, by='genes') %>% 
  gather(samplenames, norm_counts, -genes) %>% 
  ggplot() + 
  geom_point(aes(x=genes, y=log10(norm_counts), color=samplenames), position=position_jitter(w=0.1,h=0)) + 
  theme_linedraw() + 
  theme(axis.text.x=element_text(angle=90, hjust=1)) + 
  ggtitle('Plot of Log10(normalized counts) for top ten DE genes')

return(scatter_norm)
}


make_full_coldata <- function(metadata) {
sample_data <- data.frame(readr::read_csv(metadata) %>% 
               select(samplename, timepoint) %>% 
               mutate_at(vars(timepoint), factor))
return(sample_data)
}


deseq2_lrt_run <- function(counts, coldata_full, lrt_design) {

  dds_lrt <- DESeqDataSetFromMatrix(countData=counts, colData=coldata_full, design=lrt_design)
  dds_lrt$timepoint <- relevel(dds_lrt$timepoint, ref='vP0')
  dds_lrt <- DESeq(dds_lrt, test="LRT", reduced=~1)

  res_lrt <- results(dds_lrt)
  results_lrt <- c(dds_lrt_obj = dds_lrt, res_lrt_mat = res_lrt)
return(results_lrt)
}

# DESeq2 plotCounts take the unordered results, and the DDS object
plot_lrt_counts <- function(dds_lrt, res_lrt_ordered){
level_order <- c('vP0', 'vP4', 'vP7', 'vAd')
top_time <- plotCounts(dds_lrt, which.min(res_lrt_ordered$padj), intgroup=c('samplename', 'timepoint'), returnData=TRUE)
plot_lrt <- top_time %>% ggplot() + geom_point(aes(x= factor(timepoint, level=level_order), y=count, color=samplename), position=position_jitter(w=0.1,h=0)) + labs(x='Timepoint') + theme_linedraw() + ggtitle('Normalized counts over time for top DE gene from LRT analysis')
return(plot_lrt)
}

res_p7_v_p0 <- results(dds_lrt, contrast = c( "timepoint", "vP7", "vP0" ), test='Wald')
res_p7_v_p0_ordered <- res_p7_v_p0[order(res_p7_v_p0$padj),]
