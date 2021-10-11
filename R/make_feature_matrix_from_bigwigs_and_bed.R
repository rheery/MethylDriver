#' Find mean values of bigwig files for each region in a BED file
#'
#' @param bigwig_dir Path to a directory containing bigWig files.
#' @param bigwig_pattern An optional pattern used to filter the bigWig files.
#' @param genomic_regions_bed A BED file containing regions of interest
#' @param covered_bases_only Whether to only use covered bases when calculating the mean. Default is True.
#' @num_cores How many cores to use. Default is 1.
#' @return A list with the indices of the most similar n other regions for each region
make_feature_matrix_from_bigwigs_and_bed = function(bigwig_dir, bigwig_pattern = NULL, genomic_regions_bed, covered_bases_only = T, num_cores = 1){
  `%dopar%` = foreach::`%dopar%`
  `%do%` = foreach::`%do%`
  
  # Get paths to bigwig files in specified directory
  bigwig_file_paths = list.files(path = bigwig_dir, pattern = ".bw", full.names = T)
  
  # Filter for bigwig files with names matching a specified pattern
  if (!is.null(bigwig_pattern)){
    bigwig_file_paths = grep(bigwig_pattern, bigwig_file_paths, value = T)
  }
  
  # Get names of regions from BED file
  region_names = data.table::fread(genomic_regions_bed, header = F, sep = "\t", stringsAsFactors = F, select = 4)$V4
  
  # Create names for the output files that will be produced by bigWigAverageOverBed  
  average_over_bed_file_names = paste0(gsub(".bw", "", basename(bigwig_file_paths)), "_average_over_bed.txt")
  
  # Create a temporary directory to store these output files
  dir.create("temp_average_over_bed_results")
  
  # Set up cluster with specified number of cores                         
  cluster = parallel::makeCluster(num_cores, outfile = "cluster_outfile.txt")
  doParallel::registerDoParallel(cluster, cores = num_cores)
  
  # Use bigWigAverageOverBed to calculate mean value for each bigwig file for each region in the BED file and save the results in temp_average_over_bed_results directory
  foreach::foreach(bigwig_number = seq_along(bigwig_file_paths)) %dopar% {
    system2("bigWigAverageOverBed", args = c(bigwig_file_paths[bigwig_number], genomic_regions_bed, 
      paste0("temp_average_over_bed_results/", average_over_bed_file_names[bigwig_number])))
  }
  
  # Stop the cluster
  parallel::stopCluster(cluster)
  
  # Select the mean or mean0 columns from the bigWigAverageOverBed output files (corresponding to mean values for only covered bases 
  # or for all bases with uncovered bases counting as zero values, respectively) for downstream use
  mean_col = ifelse(covered_bases_only, 6, 5)  
  
  # Create a matrix with rows corresponding to genomic regions and columns corresponding to genomic features
  genomic_regions_chromatin_feature_matrix = as.matrix(foreach::foreach(result_file = 
    list.files("temp_average_over_bed_results", full.names = T), .combine = "cbind") %do% {
    data.table::fread(result_file, sep = "\t", header = F, select = mean_col)
  })
  
  # remove temp_average_over_bed_results and its contents
  unlink("temp_average_over_bed_results", recursive = T)
  
  # Assign the genomic region names as rownames of the matrix and the feature names as its column names
  rownames(genomic_regions_chromatin_feature_matrix) = region_names
  colnames(genomic_regions_chromatin_feature_matrix) = gsub("_average_over_bed.txt", "", average_over_bed_file_names)
  
  return(genomic_regions_chromatin_feature_matrix)
}
