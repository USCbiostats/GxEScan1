# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Function to convert a VCF file to a binary dosage file for GxEScan
#' 
#' Function to convert a VCF file to a binary dosage file for GxEScan. 
#' The binary formatted file is smaller that the VCF file and reads in much quicker.
#' 
#' @param vcfFilename
#' Name of VCF file
#' @param outBaseFilename
#' Base filename of output files. Three files are output.
#' A family file with extension .fam.
#' A map file with extension .bim.
#' A binary dosage file with extenstion .bdosage.
#' @param initSub
#' Amount of memory to allocate for subjects. Should be the number of subjects.
#' Additional memory will be allocated if not enough is initially allocated.
#' Having to do this can slow down the speed of the routine because of memory
#' allocations/deallocations.
#' @return
#' 0 failure
#' otherwise number of subjects read 
#' @export
VCF_to_BinaryDosage <- function(vcfFilename, outBaseFilename, initSub) {
    .Call('GxEScanR_VCF_to_BinaryDosage', PACKAGE = 'GxEScanR', vcfFilename, outBaseFilename, initSub)
}

#' Function to extract SNPs from a binary dosage file
#' 
#' Function to extract SNPs from a binary dosage file
#' 
#' @param bdosageFilename
#' Name of binary dosage file
#' @param mapFilename
#' Name of map file associated with dosage file
#' @param numSub
#' Number of subjects with data in dosage file
#' @param numSNPs
#' Number of SNPs to read in
#' @return
#' List with a vector of dosages and a matrix of probabilities
#' and a list of the input values
#' @export
ExtractDosages <- function(bdosageFilename, mapFilename, numSub, numSNPs) {
    .Call('GxEScanR_ExtractDosages', PACKAGE = 'GxEScanR', bdosageFilename, mapFilename, numSub, numSNPs)
}

#' Function to extract more SNPs from a binary dosage file
#' 
#' Function to extract more SNPs from a binary dosage file
#' 
#' @param inputs
#' List of inputs returned from ExtractDosages
#' @return
#' List with a vector of dosages and a matrix of probabilities
#' and a list of the input values
#' @export
ExtractMoreDosages <- function(inputs) {
    .Call('GxEScanR_ExtractMoreDosages', PACKAGE = 'GxEScanR', inputs)
}

#' Function to extract a SNP from a binary dosage file
#' 
#' Function to extract a SNP from a binary dosage file
#' 
#' @param bdosageFilename
#' Name of binary dosage file
#' @param mapFilename
#' Name of map file associated with dosage file
#' @param numSub
#' Number of subjects with data in dosage file
#' @param snpName
#' Name of SNP to extract
#' @param flanking
#' Number of flanking SNPs on either side to include
#' @return
#' List with a vector of dosages and a matrix of probabilities
#' and a list of the input values
#' @export
ExtractSNPDosages <- function(bdosageFilename, mapFilename, numSub, snpName, flanking) {
    .Call('GxEScanR_ExtractSNPDosages', PACKAGE = 'GxEScanR', bdosageFilename, mapFilename, numSub, snpName, flanking)
}

#' Function to merge the results from several GxEScans
#' 
#' Function to merge the results from several GxEScans for use with GxEResults
#' to produce a summary for a genome wide scan
#' 
#' @param logistic
#' Indicator if the outcome was logistic, otherwise outcome was linear
#' @param basefileNames
#' List of base file names of results
#' @param outfileName
#' Base filename for the merged results
#' @importFrom Rcpp evalCpp
#' @useDynLib GxEScanR
#' @export
GxEMerge <- function(logistic, basefileNames, outfileName) {
    .Call('GxEScanR_GxEMerge', PACKAGE = 'GxEScanR', logistic, basefileNames, outfileName)
}

#' Function to perform a GEWIS scan
#' 
#' Function to perform a GEWIS scan
#' 
#' @param y
#' Vector of outcome values, 0 or 1, all others treated as missing
#' @param x
#' Matric of covariates. Last column is tested for interaction with gene
#' @param GeneticDataFilename
#' Name of file with genetic data in a binary file (measured or dosage data)
#' @param MapFilename
#' Name of map file associated with measured genetic data file
#' @param outFilename
#' Base filename for output files
#' @param dg
#' Perform D|G test
#' @param dgxe
#' Perform Test \deqn{\beta_{x} = 0}
#' @param twodf
#' Perform 2df test betaG = 0, betaGxE = 0
#' @param threedf
#' Perform 3df test betaG = 0, betaGxE = 0, betaD = 0
#' @param ge
#' Perform G|E test
#' @param caseOnly
#' Perform G|E test on cases
#' @param controlOnly
#' Perform G|E test on cases
#' @param dgge
#' Perform 2df test betaG = 0, betaD = 0
#' @return
#' A list containint 5 matrices
#' Allele frequencies in all subjects, cases, and controls
#' The number of subjects with complete data
#' Tne number of cases with complete data
#' The parameter estimates
#' The Z statistic for each estimate
#' The 2df chi-squared statistic for beta_g = 0 and beta_GxE = 0
#' @export
ScanSNPs <- function(y, x, GeneticDataFilename, MapFilename, outFilename, dg, dgxe, twodf, threedf, ge, caseOnly, controlOnly, dgge) {
    .Call('GxEScanR_ScanSNPs', PACKAGE = 'GxEScanR', y, x, GeneticDataFilename, MapFilename, outFilename, dg, dgxe, twodf, threedf, ge, caseOnly, controlOnly, dgge)
}

#' Function to read a SNP from a measured genetic data file
#' 
#' Function to read a SNP from a measured genetic data file
#' 
#' @param BedFilename
#' Name of file with measured genetic data in plink format
#' @param numSubjects
#' Number of subjects with data in file
#' @param MapFilename
#' Name of map file associated with measured genetic data file
#' @param SNP
#' Name of SNP
#' @return
#' Chromosome Number
#' SNP Name
#' Location in base pairs
#' Dosage values
#' @export
ReadGeneticFile <- function(BedFilename, numSubjects, MapFilename, SNP) {
    .Call('GxEScanR_ReadGeneticFile', PACKAGE = 'GxEScanR', BedFilename, numSubjects, MapFilename, SNP)
}

