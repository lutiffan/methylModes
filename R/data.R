#' Illumina 450k Manifest Data
#' 
#' Annotation data for Illumina HumanMethylation450 BeadChip array (v1.2).
#' Contains probe information including chromosome, location on chromosome, 
#' position relative to CpG islands, and genetic variants under the probe.
#' Used in methylModes Shiny app to allow analysis on a subset of data.
#' Original manifest downloaded from \href{https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html}{Illumina Support}
#' 
#' @format A data.table with 485,577 rows and 6 columns:
#' \describe{
#'   \item{IlmnID}{Probe ID (e.g., "cg00000029")}
#'   \item{CHR}{Chromosome number}
#'   \item{Relation_to_UCSC_CpG_Island}{Relation to CpG island (Island, N_Shore, S_Shore, etc.)}
#'   \item{MAPINFO}{Base pair position}
#'   \item{Probe_SNPs_10}{SNPs within 10bp of probe}
#'   \item{Probe_SNPs}{SNPs 10-50bp from probe}
#' }
#' 
#' @source Illumina HumanMethylation450 BeadChip annotation
"IlluminaManifest450k"

#' Illumina EPIC Manifest Data
#' 
#' Annotation data for Illumina HumanMethylationEPIC BeadChip array (v1.0).
#' Contains probe information including chromosome, location on chromosome, 
#' position relative to CpG islands, and genetic variants under the probe.
#' Used in methylModes Shiny app to allow analysis on a subset of data.
#' Original manifest downloaded from \href{https://support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html}{Illumina Support}.
#' 
#' @format A data.table with 865,919 rows and 7 columns:
#' \describe{
#'   \item{IlmnID}{Probe ID (e.g., "cg00000029")}
#'   \item{CHR}{Chromosome number}
#'   \item{Relation_to_UCSC_CpG_Island}{Relation to CpG island (Island, N_Shore, S_Shore, etc.)}
#'   \item{MAPINFO}{Base pair position}
#'   \item{SNP_ID}{RSID of SNPs under the probe}
#'   \item{SNP_DISTANCE}{Distance from CpG to SNPs under the probe}
#'   \item{SNP_MinorAlleleFrequency}{MAF of SNPs under the probe}
#' }
#' 
#' @source Illumina HumanMethylationEPIC BeadChip annotation
"IlluminaManifestEPIC"

#' Illumina EPIC v2 Manifest Data
#' 
#' Annotation data for Illumina HumanMethylationEPIC BeadChip array (v2.0).
#' Contains probe information including chromosome, location on chromosome, 
#' position relative to CpG islands, and genetic variants under the probe.
#' Used in methylModes Shiny app to allow analysis on a subset of data.
#' SNP_DISTANCE and SNP_MinorAlleleFrequency are ommited in this copy
#' of the manifest due to file size constraints.
#' Original manifest downloaded from \href{https://support.illumina.com/downloads/infinium-methylationepic-v2-0-product-files.html}{Illumina Support}.
#' 
#' @format A data.table with 926,371 rows and 7 columns:
#' \describe{
#'   \item{IlmnID}{Probe ID (e.g., "cg25324105")}
#'   \item{CHR}{Chromosome number}
#'   \item{Relation_to_UCSC_CpG_Island}{Relation to nearest CpG island (Island, N_Shore, S_Shore, etc.)}
#'   \item{MAPINFO}{Base pair position}
#'   \item{SNP_ID}{RSID of SNPs under the probe}
#' }
#' 
#' @source Illumina HumanMethylationEPIC v2 BeadChip annotation
"IlluminaManifestEPICv2"

#' Toy Data for methylModes
#' 
#' A dataset containing randomly generated multimodal distributions. 
#' Rather than being designed to mimic DNA methylation betas distributions,
#' this data is intended to verify functionality of the Shiny app.
#' Generated using a Gaussian mixture model.
#' 
#' @format A matrix with 500 rows and 800 columns:
#' \describe{
#'   \item{rows}{Probe IDs (test1 through test500)}
#'   \item{columns}{Sample IDs (subject1 through subject800)}
#'   \item{values}{Beta values between 0 and 1}
#' }
#' 
#' @source Generated using \code{inst/generateToyMultimodalData.R}
"toyMultimodalData" 
