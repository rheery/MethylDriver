## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include = F-------------------------------------------------------
library(MethylDriver)
coad_methyldriver_results = MethylDriver:::coad_methyldriver_results

## ---- fig.align="center", echo = F--------------------------------------------
knitr::include_graphics(system.file("extdata", "methyldriver_algorithm.png", package = "MethylDriver"))

## ---- eval = F----------------------------------------------------------------
#  mean_promoter_methylation_change_coad = rowMeans(tcga_cancer_methylation_change_tables$COAD, na.rm = T)

## ---- eval = F----------------------------------------------------------------
#  promoter_neighborhoods_list = most_similar_n_regions(feature_table = epd_promoter_epigenetic_features, ranges = epd_promoters, n = 100)

## ---- eval = F----------------------------------------------------------------
#  coad_methyldriver_results = MethylDriver(region_methyl_change = mean_promoter_methylation_change_coad,
#    regions_genes = epd_promoters$gene, closest_regions_list = promoter_neighborhoods_list, multiple_correction = "fdr", q_value_cutoff = 0.05)

## -----------------------------------------------------------------------------
knitr::kable(head(coad_methyldriver_results$full_results, 10), digits = 3)

## ---- eval = F----------------------------------------------------------------
#  coad_methyldriver_results$hypermethylated_genes
#  coad_methyldriver_results$hypomethylated_genes

