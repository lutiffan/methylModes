#' Illumina 450k Manifest Data
#' 
#' Annotation data for Illumina HumanMethylation450 BeadChip array.
#' Contains probe information including chromosome, position, and CpG island relations.
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
#' Annotation data for Illumina HumanMethylationEPIC BeadChip array.
#' Contains probe information including chromosome, position, and CpG island relations.
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

#' Random Test Data for methylModes
#' 
#' A dataset containing randomly generated DNA methylation beta values
#' that mimic real methylation data. Generated using a Gaussian mixture model.
#' 
#' @format A matrix with 500 rows and 800 columns:
#' \describe{
#'   \item{rows}{Probe IDs (test1 through test500)}
#'   \item{columns}{Sample IDs (subject1 through subject800)}
#'   \item{values}{Beta values between 0 and 1}
#' }
#' 
#' @source Generated using \code{inst/generateLargeRandomData.R}
"largeRandomTestData" 