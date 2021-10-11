#' Identify similar genomic regions based on genomic features
#'
#' Uses a table with numerical features of genomic regions of interest to identify the row indices for the n most similar other regions for each region. 
#' Similarity is measured using Euclidean distance calculated using all features. A GRanges object with the genomic coordiantes of the regions can also be supplied and if so, 
#' overlapping regions will be excluded from being considered neighbours.  
#'
#' @param feature_table A data.frame or matrix with genomic features for a set of genomic regions, where features are columns and region names are row names.
#' @param ranges An optional GRanges object containing the genomic coordinates of the regions in the feature table. 
#' Row names of feature table should match names of ranges. If supplied, regions which overlap will be excluded from being considered similar to each other. 
#' @param n Number of closest neighbours to identify (default is 100).
#' @return A list with the indices of the most similar n other regions for each region
#' @export
most_similar_n_regions = function(feature_table, ranges = NULL, n = 100){
  
  # Perform a PCA on the features
  pca_result = prcomp(as.matrix(feature_table), scale. = T, center = T)$x

  # Create distance matrix from table of features
  distance_matrix = multivariance::fastdist(pca_result)
  
  # Identify overlapping regions if GRanges provided
  if(!is.null(ranges)){
    overlapping_ranges = data.frame(GenomicRanges::findOverlaps(ranges, ranges, ignore.strand = T))
    overlapping_ranges = split(overlapping_ranges$subjectHits, overlapping_ranges$queryHits)
    
    # Set distance between overlapping regions as NA so that they won't be considered as neighbours 
    for (col in 1:length(overlapping_ranges)) {
      distance_matrix[overlapping_ranges[[col]], col] = NA
    }
  }
  
  # Identify the indices of the closest n neighbour regions from the feature table for each region
  closest_regions_list = setNames(lapply(1:ncol(distance_matrix), function(x)
    order(distance_matrix[, x])[1:n]), nm = rownames(pca_result))
  
  return(closest_regions_list)
}



#' Compare methylation change of genomic regions to the change in similar regions
#'
#' #' Uses a table of numerical features for genomic regions of interest to identify the row indices for the n most similar other regions for each region. 
#' Similarity is measured using Euclidean distance calculated using all features. A GRanges object with the genomic coordiantes of the regions can be supplied and if so, 
#' overlapping regions will be excluded from being considered neighbours. 
#' 
#' @param home_region_methyl_change A named vector of methylation change values for each region, with names indicating the regions that the values refer to
#' @param closest_regions_list A list of the indices if the most similar regions for each region. Should be a list with the indices of the 
#   closest neighbour regions in covariate space for each region. Names of the list elements should be the same as those of home_region_methyl_change
#' @param multiple_correction A multiple testing correction method (default is FDR). Should be one of the choices from p.adjust.methods
#' Row names of feature table should match names of ranges. If supplied, regions which overlap will be excluded from being considered neighbours. 
#' @param regions_genes = Vector of gene names associated with each region. Can also be a list if there can be more than one gene associated with a region
#' @param q_value_cutoff = A significance threshold to use for q-values
#' @return  A list with the results of MethylDriver. Contains a vector of significant hypermethylated genes, a vector of significant hypomethylated genes and a data.frame with the name, associated gene, methylation change value, mean of the methylation change values for the most similar regions,
#' standard deviation of the methylation change values for the most similar regions, z-scores, p-values and q-values for all regions. 
#' 
#' @export
MethylDriver = function(home_region_methyl_change, regions_genes, closest_regions_list, multiple_correction = "fdr", q_value_cutoff = 0.05){
  
  # Check that names home_region_methyl_change vector has names
  if(is.null(names(home_region_methyl_change))){simpleError("home_region_methyl_change should be a named vector")}
  
  # Check that names of home_region_methyl_change match names of closest_regions_list
  if(!all(names(home_region_methyl_change) == names(closest_regions_list))){simpleError("Names of home_region_methyl_change do not match those of closest_regions_list")}
  
  # Check to see that input arguments belong to the set of allowed values
  if(!multiple_correction %in% p.adjust.methods){simpleError("Incorrect value for multiple_correction")}
  
  # Calculate for each region the distribution of the feature of interest for its neighbours
  feature_change_distribution_for_each_region = lapply(closest_regions_list, function(x) sort(home_region_methyl_change[x]))
  
  # Create a data.frame with the results. Add name of regions, genes associated with regions as a list, the value of the feature of interest for the region, 
  # the mean, the standard deviation for the feature of interest for the region's closest neighbours and the mean distance to the regions closest neighbours
  region_normalization_df = data.frame(name = names(home_region_methyl_change), gene = unname(I(regions_genes)), home_region_methyl_change = home_region_methyl_change,
    neighborhood_methyl_change_mean = sapply(feature_change_distribution_for_each_region, function(x) mean(x, na.rm = T)),
    neighborhood_methyl_change_sd = sapply(feature_change_distribution_for_each_region, sd)) 
  # Add z-score for the region by subtracting the mean of the closest neighbours from the value for a region and dividing by the standard deviation of the neighbours
  region_normalization_df$z_score = (region_normalization_df$home_region_methyl_change - region_normalization_df$neighborhood_methyl_change_mean)/region_normalization_df$neighborhood_methyl_change_sd
  # Test if the z-score differs significantly from the distribution of values for the closest neighbours (either greater than, less than or two-sided)
  region_normalization_df$p_value =  2 * pnorm(-abs(region_normalization_df$z_score))
  
  # Calculate q-value using specified method testing correction method.
  region_normalization_df$q_value = p.adjust(region_normalization_df$p_value, method = multiple_correction)
  
  # Remove closest neighbour distributions from results table 
  region_normalization_df$closest_regions_distribution = NULL
  
  # Make sure region_normalization_df is a dataframe and not a tibble
  region_normalization_df = data.frame(region_normalization_df)
  
  # Arrange results by q_value
  region_normalization_df = dplyr::arrange(region_normalization_df, q_value)
  row.names(region_normalization_df) = NULL
  
  # Get significant hypermethylated and significant hypomethylated genes
  significant_hypermethylated_genes = unique(dplyr::filter(region_normalization_df, q_value < 0.05, home_region_methyl_change > 0)$gene)
  significant_hypomethylated_genes = unique(dplyr::filter(region_normalization_df, q_value < 0.05, home_region_methyl_change < 0)$gene)
  
  # Return specified results object
  return(list(hypermethylated_genes = significant_hypermethylated_genes, hypomethylated_genes = significant_hypomethylated_genes, full_results = region_normalization_df))
}
