#' Promoters of autosomal protein-coding genes
#'
#' A GRanges object containing promoters of autosomal protein-coding genes  from Eukaryotic Promoter Database (EPD). 
#' The 1 KB upstream of the transcription start sites were used to define the promoters.
#'
#' @format A GRanges object with 23,208 ranges with metadata columns giving the name of the promoter and the name of the associated gene
#' @source Promoters were downloaded from \link{https://epd.epfl.ch/}
"epd_promoters"

#' Epigenetic features at EPD promoters from Roadmap Epigenomics Project
#'
#' A data.frame with epigenetic features for \code{epd_autosomal_protein_coding_gene_promoters} from the Roadmap Epigenomics project
#'
#' @format A data.frame with 23,208 rows and 878 columns, with rows corresponding to promoters and columns to epigenetic features from the Roadmap Epigenetics Project
#' @source Roadmap data was downloaded from \link{https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/} and promoter values were calculated using the 
#' bigwigaverageoverbed tool from UCSC
"epd_promoter_epigenetic_features"

#' Most similar 100 promoters for each promoter
#'
#' A list with the indices of the 100 most similar promoters for each promoter in \code{epd_autosomal_protein_coding_gene_promoters} to be used with \code{MethylDriver}
#'
#' @format A named list of length 23,208 where the names of the list elements correspond to the promoter names from \code{epd_autosomal_protein_coding_gene_promoters} and the 
#' values to the indices of the 100 most similar promoters in \code{epd_autosomal_protein_coding_gene_promoters} for that promoter. 
"promoter_neighborhoods_100"

#' Mean methylation change at EPD promoters in TCGA cancer types
#'
#' A list with the methylation change at each promoter from \code{epd_promoters} in matching normal-tumour samples for 13 TCGA cancer types:
#' BLCA = Bladder Urothelial Carcinoma, BRCA = Bladder Urothelial Carcinoma, COAD = Colon Adenocarcinoma, ESCA = Esophageal Carcinoma, HNSC = Head and Neck Squamous Cell Carcinoma, 
#' KIRC = Kidney Renal Clear Cell Carcinoma, KIRP = Kidney renal papillary cell carcinoma, LIHC = Hepatocellular Carcinoma, LUAD = Lung Adenocarcinoma, 
#' LUSC = Lung Squamous Cell Carcninoma, PRAD = Prostate Adenocarcinoma, THCA = Thyroid Cancer, UCEC = Uterine Corpus Endometrial Carcinoma
#'
#' @format A list of length 13 where each element is a table 23,208 rows and a column with the methylation change for each matching normal-tumour pair
"tcga_cancer_methylation_change_tables"