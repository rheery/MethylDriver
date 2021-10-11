# Create a matrix with epigenetic features for promoters and create a list with promoter neighborhoods of 100

# Load required packages
library(MethylDriver)
library(BSgenome.Hsapiens.UCSC.hg19)
library(foreach)
library(doParallel)
cluster = makeCluster(5, outfile = "cluster_outfile.txt")
registerDoParallel(cluster)

# All files in system.file("extdata", "all_roadmap_bigwigs.txt", package = "MethylDriver") should be downloaded and placed in a directory called roadmap_bigwigs

# Next a bigWig with the position of CpG sites in the genome is created

# Create empty RLE for lengths of all chromosomes
chrom_lengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[1:23]
chrom_ranges = setNames(GRanges(seqnames = names(chrom_lengths), ranges = IRanges(start = 1, end = chrom_lengths)), nm = names(chrom_lengths))
chrom_rles = mclapply(chrom_lengths, function(x) Rle(rep(0, x)), mc.cores = 5)
chrom_rles = RleList(chrom_rles, compress = F)

# Create views object with location of all CpGs
params  = new("BSParams", X = BSgenome.Hsapiens.UCSC.hg19, FUN= matchPattern, exclude = c("_"))
all_cpgs  = bsapply(BSParams = params, pattern = "CG")  

# Make GRanges object from views object
cpg_genome_ranges_hg19 = GRanges(seqnames = rep(names(all_cpgs),  lengths(all_cpgs)), ranges = IRanges(unlist(sapply(all_cpgs, start)), unlist(sapply(all_cpgs, start))))

# Remove objects no longer needed
rm(params, all_cpgs)

# Make a dataframe from the GRanges
cpg_genome_ranges_hg19_df = data.frame(cpg_genome_ranges_hg19, stringsAsFactors = F)

# Convert dataframe into RLE
chrom_rles = setNames(foreach (chrom = names(chrom_rles), .packages = c("dplyr", "BSgenome.Hsapiens.UCSC.hg19"))  %dopar% {
  chrom_table = dplyr::filter(cpg_genome_ranges_hg19_df, seqnames == chrom)
  chrom_rles[[chrom]][chrom_table$start] = 1
  chrom_rles[[chrom]]
}, nm = names(chrom_rles))

chrom_rles = RleList(chrom_rles, compress = F)

# Save CpG chrom RLE and export as a BigWig
export.bw(chrom_rles, "roadmap_bigwigs/hg19_cpgs.bw")

# The appropriate version of bigWigAverageOverBed should be downloaded from https://hgdownload.soe.ucsc.edu/admin/exe/ and made available on the search path.

# Create PCA feature value matrix for promoter ranges autosomal overlapping probes. May take a while so raise num_cores to increase running speed. 
epd_promoter_feature_matrix = MethylDriver:::make_feature_matrix_and_pca_from_bigwigs_and_bed(bigwig_dir = "roadmap_bigwigs/", 
  genomic_regions_bed = system.file("extdata", "epd_promoters.bed", package = "MethylDriver"), covered_bases_only = T, num_cores = 4)

# Create a list with promoter neighborhoods of size 100
promoter_neighborhoods_100 = most_similar_n_regions(feature_table = epd_promoter_feature_matrix, ranges = epd_promoters, n = 100)

saveRDS(promoter_neighborhoods_100, "promoter_neighborhoods_100.rds")