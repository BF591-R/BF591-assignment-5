library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param subset(list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, subset) {
  meta <- read_csv(metafile_csv) %>% 
    dplyr::select(samplename, timepoint) %>%
    filter(timepoint %in% subset) %>%
    mutate(timepoint = factor(timepoint, levels=subset))
  
  subset_samples <- pull(meta, samplename)
  counts <- data.matrix(read.delim(counts_csv, sep='\t', row.names='gene'))
  subset <- counts[, subset_samples]
  
  se <- SummarizedExperiment(assays = list(counts = subset),
                             colData = meta)
  metadata(se)$model <- counts ~ timepoint
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  res <- results(dds) %>% as.data.frame()
  
  return(list('res' = res, 'dds' = dds))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  labeled <- deseq2_res %>%
    as_tibble(rownames='genes') %>%
    mutate(volc_plot_status = case_when(log2FoldChange > 0 & padj < padj_threshold ~ 'UP', 
                                        log2FoldChange < 0 & padj < padj_threshold ~ 'DOWN', 
                                        TRUE ~ 'NS'))
  return(labeled)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  pval_plot <- labeled_results %>% 
    ggplot(aes(pvalue)) + 
    geom_histogram(bins=50, color='black', fill='lightblue') + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle('Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)')
  
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  logfc_plot <- labeled_results %>% 
    filter(padj < padj_threshold) %>% 
    ggplot(aes(log2FoldChange)) + 
    geom_histogram(bins=100, color='black', fill='light blue') + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle('Histogram of Log2FoldChanges for DE Genes (vP0 vs. vAd)')
  return(logfc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  top_genes <- labeled_results %>% 
    slice_min(padj, n = num_genes) %>% 
    dplyr::select(genes)
  
  dds_obj <- estimateSizeFactors(dds_obj)
  norm_counts <- counts(dds_obj, normalized=TRUE) %>% as_tibble(rownames='genes')
  
  scatter_norm <- top_genes %>% 
    left_join(norm_counts, by='genes') %>% 
    gather(samplenames, norm_counts, -genes) %>% 
    ggplot() + 
    geom_point(aes(x=genes, y=log10(norm_counts), color=samplenames), position=position_jitter(w=0.1,h=0)) + 
    theme_minimal() +  
    theme(axis.text.x=element_text(angle=90, hjust=1)) + 
    xlab(' ') +
    ggtitle('Plot of Log10(normalized counts) for top ten DE genes')
  
  return(scatter_norm)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  volcano_plot <- labeled_results %>% 
    ggplot() + 
    geom_point(mapping=aes(x=log2FoldChange, y=-log10(padj), color=volc_plot_status)) + 
    geom_hline(yintercept = -log10(0.1), linetype = "dashed")  + 
    theme_minimal() +
    ggtitle('Volcano plot of DESeq2 differential expression results (vP0 vs. vAd)')
  return(volcano_plot)
}

#' Function to run fgsea on DESeq2 results
#'
#' @param labeled_results (tibble): the labeled results from DESeq2
#' @param gmt (str): the path to the GMT file
#' @param min_size: the threshold for minimum size of the gene set
#' @param max_size: the threshold for maximum size of the gene set
#'
#' @return tibble containing the results from running fgsea using descending
#' log2foldchange as a ranking metric
#' @export
#'
#' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
run_gsea <- function(labeled_results, gmt, min_size, max_size) {
  
  labeled_results <- labeled_results %>% separate(genes, sep='\\.', into='genes', remove=TRUE)
  gene_ids <- labeled_results %>% pull(genes)
  
  human <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
  mouse <- useMart('ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')
  
  hgnc_symbols <- getLDS(attributes=c('ensembl_gene_id'), 
                         filters='ensembl_gene_id', 
                         values=gene_ids, 
                         mart = mouse, 
                         attributesL = c('hgnc_symbol'), 
                         martL = human, 
                         uniqueRows= TRUE)
  
  hgnc_results <- labeled_results %>% left_join(hgnc_symbols, by=c('genes' = 'Gene.stable.ID'))
  
  rnks <- hgnc_results %>% 
    drop_na(HGNC.symbol, log2FoldChange) %>% 
    distinct(HGNC.symbol, log2FoldChange, .keep_all=TRUE) %>%
    arrange(desc(log2FoldChange)) %>% 
    dplyr::select(HGNC.symbol, log2FoldChange) %>% 
    deframe()
  
  c2_pathways <- gmtPathways('c2.cp.v7.5.1.symbols.gmt')
  
  fgsea_results <- fgsea(c2_pathways, rnks, minSize=min_size, maxSize=max_size) %>% as_tibble()
  
  return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  
  top_pos <- fgsea_results %>% slice_max(NES, n=num_paths) %>% pull(pathway)
  top_neg <- fgsea_results %>% slice_min(NES, n=num_paths) %>% pull(pathway)
  
  subset <- fgsea_results %>% 
    filter(pathway %in% c(top_pos, top_neg)) %>%
    mutate(pathway = factor(pathway)) %>%
    mutate(plot_name = str_replace_all(pathway, '_', ' '))
  
  plot <- subset %>% 
    mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
    ggplot() +
    geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
    scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
    theme_minimal(base_size = 8) +
    ggtitle('fgsea results for Hallmark MSigDB gene sets') +
    ylab('Normalized Enrichment Score (NES)') +
    xlab('') +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
    coord_flip()
  return(plot)
}
