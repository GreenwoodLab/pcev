#' pcev: A package for computing principal components of explained variance.
#'
#' PCEV is a statistical tool for the analysis of a mutivariate response vector.
#' It is a dimension-reduction technique, similar to Principal Components
#' Analysis (PCA), which seeks the maximize the proportion of variance (in the
#' response vector) being explained by a set of covariates.
#'
#' @section pcev functions:
#' \code{\link{estimatePcev}}
#' \code{\link{computePCEV}}
#' \code{\link{PcevObj}}
#' \code{\link{permutePval}}
#' \code{\link{wilksPval}}
#' \code{\link{roysPval}}
#'
#' @docType package
#' @name pcev-package
NULL

####################################
# Documentation for the datasets----

#' Methylation values around BLK gene
#' 
#' A dataset containing methylation values for cell-separated samples. The methylation was measured
#' using bisulfite sequencing. The data also contains the genomic position of these CpG sites, as
#' well as a binary phenotype (i.e. whether the sample comes from a B cell).
#' 
#' Methylation was first measured at 24,068 sites, on 40 samples. Filtering was performed to keep
#' the 25\% most variable sites. See the vignette for more detail.
#' 
#' A second sample of the methylation dataset was extracted. This second dataset contains
#' methylation values at 1000 CpG dinucleotides.
#' 
#' @format The data comes in four objects: 
#' \describe{ 
#' \item{methylation}{Matrix of methylation
#'   values at 5,986 sites measured on 40 samples} 
#' \item{pheno}{Vector of phenotype, indicating
#'   whether the sample comes from a B cell} 
#' \item{position}{Data frame recording the position of
#'   each CpG site along the chromosome} 
#' \item{index}{Index vector used in the computation of
#'   PCEV-block} 
#' \item{methylation2}{Matrix of methylation values at 1000 sites measured on 40
#'   samples} 
#' \item{pheno2}{Vector of phenotype, indicating the cell type of the sample (B cell, T cell, or Monocyte)} 
#' \item{position2}{Data frame recording the position of each CpG site along the chromosome} }
#' @source Tomi Pastinen, McGill University, GEnome Quebec.
"methylation"

#' @rdname methylation
"pheno"

#' @rdname methylation
"position"

#' @rdname methylation
"index"

#' @rdname methylation
"pheno2"

#' @rdname methylation
"position2"

#' @rdname methylation
"methylation2"
