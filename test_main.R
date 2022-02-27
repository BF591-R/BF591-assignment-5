source('main.R')
library('testthat')


test_names <- paste0('ENSMUSG', 0:5)

test_data <- data.frame(
  log2FoldChange = c(rep(2, 3), rep(-2, 3)),
  pval = rep(c(.001, .01, .1), 2),
  padj = rep(c(.01, .10, .15), 2))

row.names(test_data) <- test_names

test_tib <- test_data %>% 
  as_tibble(rownames='genes') %>%
  mutate(volc_plot_status = c('UP', rep('NS', 2), 'DOWN', rep('NS', 2)))


test_that('SummarizeExperiments object holds the right data', {
  se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
  
  expect_true(class(se) == 'SummarizedExperiment')
  expect_equal(dim(assays(se)$counts), c(55416, 4))
  expect_setequal(colnames(se), c('vP0_1', 'vP0_2', 'vAd_1', 'vAd_2'))
  expect_true(is.matrix(assays(se)$counts), 'DESeq2 expects a matrix of counts as an input')
  
})

test_that('return_deseq_res returns a list of results and the dds object', {
  
  se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
  deseq2_results <- return_deseq_res(se, ~timepoint)
  
  expect_equal(length(deseq2_results), 2)
  expect_equal(class(deseq2_results), 'list')
  
})

test_that('label_res correctly labels test data with known volc_plot_status', {
  
  test_subset <- test_tib %>% dplyr::select(log2FoldChange, padj, volc_plot_status)
  

  func_labels <- label_res(test_data, .10)
  func_cols <- func_labels %>% dplyr::select(log2FoldChange, padj, volc_plot_status)
  func_cols
  
  expect_equal(func_cols, test_subset)
  
})

test_that('DESeq2 results return the same directionality as the most significant gene', {
  se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
  dds <- DESeqDataSet(se, design = ~timepoint)
  dds <- DESeq(dds)
  res <- results(dds) %>% as_tibble(rownames='genes')
  
  fc_sign_neg <- res %>% filter(genes == 'ENSMUSG00000026418.17') %>% dplyr::select(log2FoldChange) %>% pull()
  
  fc_sign_pos <- res %>% filter(genes == 'ENSMUSG00000002500.16') %>% dplyr::select(log2FoldChange) %>% pull()
  
  
  expect_true(fc_sign_neg < 0, 'Check the reference level in DESeq2')
  expect_true(fc_sign_pos > 0, 'Check the reference level in DESeq2')
})

test_that("plot pvals and log2fc return histograms", {


  
  pval_plot <- plot_pvals(test_tib)
  log2fc_plot <- plot_log2fc(test_tib, .10)
  
  # testing all geoms in the ggplot object
  geoms_pval <- c()
  for (geom in pval_plot$layers) {
    geoms_pval <- c(geoms_pval, class(geom$geom)[1])
  }
  
  geoms_log2fc <- c()
  for (geom in log2fc_plot$layers) {
    geoms_log2fc <- c(geoms_log2fc, class(geom$geom)[1])
  }
  
  expect_true("GeomBar" %in% geoms_pval)
  expect_true("GeomBar" %in% geoms_log2fc)
})

test_that('volcano plot has the right labels', {
  volc_plot <- plot_volcano(test_tib)

  expect_true(volc_plot$labels$colour == 'volc_plot_status')
  
  geoms_volc <- c()
  for (geom in volc_plot$layers) {
    geoms_volc <- c(geoms_volc, class(geom$geom)[1])
  }
  
  expect_true("GeomPoint" %in% geoms_volc)
})

test_that('fgsea results have the right results columns and is tibble', {

  se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
  dds <- DESeqDataSet(se, design = ~timepoint)
  dds <- DESeq(dds)
  res <- results(dds) %>% as_tibble(rownames='genes')
  
  fgsea_res <- run_gsea(res, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
  
  fgsea_plot <- top_pathways(fgsea_res, 10)
  
  expect_true(sum(fgsea_plot$data$NES > 0) == 10)
  expect_true(sum(fgsea_plot$data$NES > 0) == 10)
  
  expect_true(is_tibble(fgsea_res))
  expect_true(all(c('NES', 'ES', 'pval', 'padj', 'log2err', 'size', 'leadingEdge') %in% names(fgsea_res)))
  
})

