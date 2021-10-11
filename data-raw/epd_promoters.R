# Process the BED file from EPDnew for hg19 to make promoters GRanges objects for protein-coding genes which overlap Infinium 450k probes

# Load required packages
library(dplyr)
library(rtracklayer)
library(HGNChelper)

# Download HGNC complete table. Note that the table from 2020-09-01 was used for our analysis so there may be some very slihgt differences
hgnc_full_table = read.csv("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/quarterly/tsv/hgnc_complete_set_2020-10-01.txt", sep = "\t", header = T)

# Create tables for protein-coding genes
protein_coding_genes = filter(hgnc_full_table, locus_type == "gene with protein product", status == "Approved")$symbol

# Download EPDnew Version 6 hg19 promoters 
epd_promoters = import.bed("ftp://ccg.epfl.ch/epdnew/H_sapiens/006/Hs_EPDnew_006_hg19.bed", genome = "hg19")

# Sort promoter granges
epd_promoters = sort(epd_promoters, ignore.strand = T)

# Promoters are resized so that they are 1000 bp by adding necessary number of bases upstream
epd_promoters = resize(epd_promoters, width = 1000, fix = "end")

# Score and thick columns are removed since they are not needed
epd_promoters$score = NULL
epd_promoters$thick = NULL

# Extract gene names from promoter names
epd_promoters$gene = gsub("_.", "", epd_promoters$name)

# Usew HGNChelper to check that gene symbols are up to date. Note that the current map on the 1/9/20 was downloaded using getCurrentHumanMap() and used as the map argument for checkGeneSymbols. 
# Using a different map may give different gene names for a small number of genes. 
epd_promoters$gene = checkGeneSymbols(epd_promoters$gene)$Suggested.Symbol

# Genes which could not be associated with a HGNC symbol or matched multiple symbols or were annotated on the wrong chromosome (according to HGNC) are removed
epd_promoters$gene[grep("///", epd_promoters$gene)] = NA
epd_promoters = epd_promoters[!is.na(epd_promoters$gene)]

# Subset for genes annotated as protein-coding by HGNC. 
epd_pcg_promoters =  epd_promoters[epd_promoters$gene %in% protein_coding_genes]

# Subset for autosomal protein-coding genes. 
epd_pcg_promoters_autosomal = epd_pcg_promoters[seqnames(epd_pcg_promoters) %in% paste0("chr", 1:22)]

# Read in Infinium hg19 probe GRanges which correspond to good probes from GDC 450k tables
illumina_probe_hg19_ranges = import.bed(system.file("extdata", "illumina_450k_probe_granges.bed", package = "MethylDriver"))

# Find promoters which overlap probes from the Infinium 450K array. 
epd_pcg_promoters_autosomal_overlapping_probes = subsetByOverlaps(epd_pcg_promoters_autosomal, illumina_probe_hg19_ranges) 

# Save epd_pcg_promoters_autosomal_overlapping_probes as an RDS file
saveRDS(epd_pcg_promoters_autosomal_overlapping_probes, "epd_promoters.rds")
export.bed(epd_pcg_promoters_autosomal_overlapping_probes, "epd_promoters.bed")