#' @importFrom Rcpp evalCpp
#' @importFrom data.table fread setkey data.table
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis lines plot points rect title
#' @importFrom stats qqplot
#' @importFrom utils read.table
#' @useDynLib GxEScanR
NULL

#' Function to perform a GEWIS scan
#' 
#' Function to perform a GEWIS scan
#' 
#' @param familyFile
#' Name of file with pedigree data in plink format
#' @param continuous
#' Indicator if outcome is continous
#' @param zeroOne
#' Indicator if phenotype is coded 0,1. Default is 1,2.
#' @param pmissing
#' Value to be treated as missing for phenotype
#' @param covariateFile
#' Name of file containing covariate data in plink format
#' @param header
#' Indicator if there is a header in the covariate file
#' @param covsUsed
#' List of covaraites to be used
#' @param interactionCov
#' Covariate to be used in interaction
#' @param cmissing
#' Value of covariate to be considered missing
#' @param useSex
#' Indicator if value for sex from family file is to be used as a covariate
#' @param sexInteraction
#' Indicator if value for sex from family file is to be used as the interacting covariate
#' @param phenotype
#' Value from covariate file to be used as phenotype. If the value is NA, the value from
#' the family file will be used as the phenotype
#' @param geneticFile
#' Name of binary genetic data file may be measured or dosages
#' @param mapFile
#' Name of map file associated with genetic data file
#' @param outfile
#' Base filename for output files
#' @param dg
#' Perform gene only scan
#' @param dgxe
#' Perform scan with gene by environment interaction in model
#' @param twodf
#' Perform two degree of freedom test (g and gxe)
#' @param threedf
#' Perform three degree of freedom test
#' @param ge
#' Perform scan of g given e
#' @param caseOnly
#' Fit case only model
#' @param controlOnly
#' Fit control only model
#' @param dgge
#' Perfromf two degree of freedom test (d|g, g|e)
#' @export 
GxEScan <- function(familyFile, continuous = FALSE, zeroOne = FALSE, pmissing = NA,
                    covariateFile, header = TRUE, covsUsed, interactionCov, cmissing = NA,
                    useSex = FALSE, sexInteraction = FALSE, phenotype = NA,
                    geneticFile, mapFile, outfile = "",
                    dg = TRUE, dgxe = TRUE, twodf = TRUE, threedf = TRUE,
                    ge = TRUE, caseOnly = TRUE, controlOnly = TRUE, dgge = TRUE) {
  fam <- ReadFamily(filename = familyFile, continuous = continuous, zeroOne = zeroOne, pmissing = pmissing)
  if (missing(covariateFile) == FALSE) { 
    if (missing(covsUsed)) {
      if (missing(interactionCov))
        cov <- ReadCovariates(filename = covariateFile, header = header, cmissing = cmissing)
      else
        cov <- ReadCovariates(filename = covariateFile, header = header, interaction = interactionCov, cmissing = cmissing)
    } else {
      if (missing(interactionCov))
        cov <- ReadCovariates(filename = covariateFile, header = header, used = covsUsed, cmissing = cmissing)
      else
        cov <- ReadCovariates(filename = covariateFile, header = header, used = covsUsed, interaction = interactionCov, cmissing = cmissing)
    }
    gxedata <- CreateGxEScanDataset(family = fam, covariates = cov, sex = useSex, sexInteraction = sexInteraction, phenotype = phenotype)
  } else {
    gxedata <- CreateGxEScanDataset(family = fam, sex = useSex, sexInteraction = sexInteraction, phenotype = phenotype)
  }
  
  return (ScanSNPs(y = gxedata$Y, x = gxedata$X, GeneticDataFilename = geneticFile, MapFilename = mapFile, outFilename = outfile,
                   dg = dg, dgxe = dgxe, twodf = twodf, threedf = threedf, ge = ge, caseOnly = caseOnly,
                   controlOnly = controlOnly, dgge = dgge))
}